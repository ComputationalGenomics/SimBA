// ============================================================================
// SimBA-hap
//
// Copyright (c) 2016 IBM Research
// ============================================================================

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>

#include <boost/multi_array.hpp>

#include <lemon/lp.h>

// ============================================================================
// View
// ============================================================================

// ----------------------------------------------------------------------------
// Class reference_view
// ----------------------------------------------------------------------------

template <class T>
class reference_view
{
    T* value_;

public:
    reference_view() :
        value_()
    {}

    reference_view(T & value) :
        value_(std::addressof(value))
    {}

    T & get() const
    {
        SEQAN_ASSERT(value_);
        return *value_;
    }

    operator T & () const
    {
        SEQAN_ASSERT(value_);
        return *value_;
    }

    reference_view<T> operator=(T & value)
    {
        value_ = std::addressof(value);
        return *this;
    }
};

// ----------------------------------------------------------------------------
// Function view
// ----------------------------------------------------------------------------

template <class T>
inline reference_view<T> view(T & t)
{
    return reference_view<T>(t);
}

// ============================================================================
// Extensions to STL
// ============================================================================

// ----------------------------------------------------------------------------
// Operator <<
// ----------------------------------------------------------------------------

namespace std
{
template <typename T1, typename T2>
std::ostream& operator<<(std::ostream & s, std::pair<T1, T2> const & p)
{
    s << '<' << std::get<0>(p) << ',' << std::get<1>(p) << '>';
    return s;
}

template <typename T>
std::ostream& operator<<(std::ostream & s, std::vector<T> const & v)
{
    s << "[ ";
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(s, " "));
    s << "]";
    return s;
}

template <typename T, size_t SIZE>
std::ostream& operator<<(std::ostream & s, std::array<T, SIZE> const & v)
{
    s << "[ ";
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(s, " "));
    s << "]";
    return s;
}
}

namespace boost {
template <typename T, size_t SIZE>
std::ostream& operator<<(std::ostream & s, boost::multi_array<T, SIZE> const & m)
{
    s << "[ ";
    std::copy(m.data(), m.data() + m.num_elements(), std::ostream_iterator<T>(s, " "));
    s << "]";
    return s;
}
}

// ============================================================================
// VCF I/O
// ============================================================================

// ----------------------------------------------------------------------------
// Function read_alleles()
// ----------------------------------------------------------------------------
// Read ref/alt strings.

template <typename alleles_t>
inline void read_alleles(alleles_t & alleles, seqan::VcfRecord const & record)
{
//    front(alleles) = record.ref;
//    clear(back(alleles));
//    auto alleles_it = directionIterator(record.alt, seqan::Input());
//    seqan::readUntil(back(alleles), alleles_it, seqan::EqualsChar<','>());

    clear(alleles);
    appendValue(alleles, record.ref);
    strSplit(alleles, record.alt, seqan::EqualsChar<','>());
}

// ----------------------------------------------------------------------------
// Function read_genotype()
// ----------------------------------------------------------------------------
// Read genotype string.

template <typename genotype_t, typename genotype_info_t>
inline void read_genotype(genotype_t & genotype, genotype_info_t const & genotype_info)
{
    genotype.clear();

    // Read first field, assuming it is GT.
    auto genotype_it = directionIterator(genotype_info, seqan::Input());
    seqan::readUntil(genotype, genotype_it, seqan::EqualsChar<':'>());

    // Convert genotype string to vector by removing slashes.
    genotype.erase(std::remove(genotype.begin(), genotype.end(), '/'), genotype.end());
    genotype.erase(std::remove(genotype.begin(), genotype.end(), '|'), genotype.end());
}

// ----------------------------------------------------------------------------
// Function write_genotype()
// ----------------------------------------------------------------------------
// Write genotype string.

template <typename genotype_info_t, typename genotype_t>
inline void write_genotype(genotype_info_t & genotype_info, genotype_t const & genotype)
{
    seqan::clear(genotype_info);

    std::for_each(genotype.begin(), genotype.end() - 1, [&genotype_info](auto allele)
    {
        seqan::appendNumber(genotype_info, allele.get());
        seqan::appendValue(genotype_info, '|');
    });

    seqan::appendNumber(genotype_info, genotype.back().get());
}

// ============================================================================
// Genotypes
// ============================================================================

// ----------------------------------------------------------------------------
// Function get_ploidy()
// ----------------------------------------------------------------------------

template <typename genotype_t>
inline uint32_t get_ploidy(genotype_t const & genotype)
{
    return genotype.size();
}

// ----------------------------------------------------------------------------
// Function is_unknown()
// ----------------------------------------------------------------------------

template <typename genotype_t>
inline bool is_unknown(genotype_t const & genotype)
{
    return genotype.back() == '.';
}

// ----------------------------------------------------------------------------
// Function alt_dosage()
// ----------------------------------------------------------------------------

//template <typename alt_t>
//struct alt_allele
//{
//    static const bool value = true;
//};
//
//template <>
//struct alt_allele<char>
//{
//    static const char value = '1';
//};

template <typename genotype_t>
inline uint32_t alt_dosage(genotype_t const & genotype)
{
    return std::count(genotype.begin(), genotype.end(), '1');
//    return std::count(genotype.begin(), genotype.end(), alt_allele<typename genotype_t::value_type>::value);
}

// ============================================================================
// Dosages
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction dosages<>::type
// ----------------------------------------------------------------------------

template <typename value_t, uint32_t n_ploidy>
struct dosages
{
    typedef std::array<value_t, n_ploidy + 1u> type;
};

// ----------------------------------------------------------------------------
// Function alt_dosages()
// ----------------------------------------------------------------------------

template <typename value_t, uint32_t n_ploidy, typename genotypes_t>
inline typename dosages<value_t, n_ploidy>::type alt_dosages(genotypes_t const & genotypes)
{
    typename dosages<value_t, n_ploidy>::type dosages;

    std::fill(dosages.begin(), dosages.end(), 0u);

    for (auto const & genotype : genotypes)
        dosages[alt_dosage(genotype)]++;

    return dosages;
}

// ============================================================================
// Markers
// ============================================================================

// ----------------------------------------------------------------------------
// Class markers
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
struct markers
{
    typedef std::pair<uint32_t, uint32_t> position_t;
//    typedef std::array<seqan::IupacString, 2> alleles_t;
    typedef seqan::StringSet<seqan::IupacString, seqan::Owner<seqan::ConcatDirect<>>> alleles_t;
    typedef typename dosages<float, n_ploidy>::type dosages_t;

    std::vector<position_t> positions; // positions[m] = (chr, pos)
    std::vector<alleles_t> alleles;    // alleles[m] = (ref, alt)
    std::vector<dosages_t> dosages;    // dosages[m] = [f(0),f(1),...,f(p)]

    inline void read(std::string const & vcf_filename_in);

    template <typename generator_t>
    inline void simulate(uint32_t n_samples, uint32_t n_markers, generator_t && generator);

    inline void normalize(uint32_t n_samples);
};

// ----------------------------------------------------------------------------
// Method markers::read()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
inline void markers<n_ploidy>::read(std::string const & vcf_filename_in)
{
    typedef markers::position_t position_t;
    typedef markers::alleles_t alleles_t;
    typedef markers::dosages_t dosages_t;
    typedef std::vector<char> genotype_t;

    position_t position;
    alleles_t alleles;
    dosages_t dosages;
    genotype_t genotype;

    seqan::VcfHeader vcf_header;
    seqan::VcfRecord vcf_record;

    seqan::VcfFileIn vcf_file_in(vcf_filename_in.c_str());
    readHeader(vcf_header, vcf_file_in);

    while (!atEnd(vcf_file_in))
    {
        readRecord(vcf_record, vcf_file_in);

        position_t position = std::make_pair(vcf_record.rID, vcf_record.beginPos);
        read_alleles(alleles, vcf_record);

        if (length(alleles) > 2)
        {
            std::cerr << "INPUT VARIANT @ " << position << " POLYALLELIC" << std::endl;
            continue;
        }

//        dosages.resize(n_ploidy);
        std::fill(dosages.begin(), dosages.end(), 0u);

        // Read all sample genotypes.
        for (auto const & genotype_info : vcf_record.genotypeInfos)
        {
            read_genotype(genotype, genotype_info);

            // Skip unknown genotypes.
            if (is_unknown(genotype))
            {
                std::cerr << "INPUT GENOTYPE @ " << position << " UNKNOWN" << std::endl;
                continue;
            }

            // Stop on wrong ploidy.
            if (get_ploidy(genotype) != n_ploidy)
                throw seqan::IOError("Input ploidy does not match VCF genotypes");

            dosages[alt_dosage(genotype)]++;
        }

        this->positions.push_back(position);
        this->alleles.push_back(alleles);
        this->dosages.push_back(dosages);

        std::cerr << "INPUT DOSAGES @ " << position << " # " << dosages << std::endl;
    }

    std::cerr << "=================================================================" << std::endl << std::endl;
}

// ----------------------------------------------------------------------------
// Method markers::simulate()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy> template <typename generator_t>
inline void markers<n_ploidy>::simulate(uint32_t n_samples, uint32_t n_markers, generator_t && generator)
{
    seqan::ignoreUnusedVariableWarning(n_samples);
    seqan::ignoreUnusedVariableWarning(n_markers);
    seqan::ignoreUnusedVariableWarning(generator);
}

// ----------------------------------------------------------------------------
// Method markers::normalize()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
inline void markers<n_ploidy>::normalize(uint32_t n_samples)
{
    std::for_each(dosages.begin(), dosages.end(), [this, n_samples](auto & marker_dosages)
    {
        auto dosages_sum = std::accumulate(marker_dosages.begin(), marker_dosages.end(), 0.0f);
        std::transform(marker_dosages.begin(), marker_dosages.end(), marker_dosages.begin(),
                      [n_samples, dosages_sum](auto d) { return n_samples * d / dosages_sum; });
    });
}

// ============================================================================
// Population
// ============================================================================

// ----------------------------------------------------------------------------
// Class population
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
struct population
{
private:
    typedef std::vector<uint16_t> founder_alt_t;
    typedef std::array<reference_view<uint16_t>, n_ploidy> sample_alt_t;
    typedef std::vector<sample_alt_t> samples_alt_t;

public:
    typedef std::vector<uint32_t> founders_freqs_t;         // [founder]
    typedef boost::multi_array<uint32_t, 2> founders_map_t; // [sample][haplotype]
    typedef std::vector<founder_alt_t> founders_alts_t;     // [marker][founder]
    typedef std::vector<samples_alt_t> samples_alts_t;      // [marker][sample][haplotype]

    inline void resize(uint32_t n_founders, uint32_t n_samples, uint32_t n_markers);

    template <typename generator_t>
    inline void simulate(generator_t && generator);

    inline void write(std::string const & vcf_filename_out);

    population(markers<n_ploidy> const & markers_in) :
        markers_in(markers_in),
        n_founders(),
        n_samples(),
        n_markers()
    {}

private:
    markers<n_ploidy> const & markers_in;

    founders_freqs_t founders_freqs;
    founders_map_t founders_map;
    founders_alts_t founders_alts;
    samples_alts_t samples_alts;

    uint32_t n_founders;
    uint32_t n_samples;
    uint32_t n_markers;

    template <typename generator_t>
    inline void simulate_founders_freqs(generator_t && generator);

    template <typename generator_t>
    inline void simulate_founders_map(generator_t && generator);

    template <typename generator_t>
    inline void simulate_founders_alts(generator_t && generator);

    inline void fit_founders_alts();
    inline void fit_founders_alts(uint32_t marker_id);

    inline void assign_samples_alts();
};

// ----------------------------------------------------------------------------
// Method population::resize()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
inline void population<n_ploidy>::resize(uint32_t n_founders, uint32_t n_samples, uint32_t n_markers)
{
    SEQAN_ASSERT_GT(n_samples, 0u);
    SEQAN_ASSERT_LEQ(n_founders, n_samples * n_ploidy);
    SEQAN_ASSERT_GT(n_markers, 0u);

    this->n_founders = n_founders;
    this->n_samples = n_samples;
    this->n_markers = n_markers;

    founders_freqs.resize(n_founders, 1u);
    founders_map.resize(boost::extents[n_samples][n_ploidy]);
    founders_alts.resize(n_markers, founder_alt_t(n_founders));
    samples_alts.resize(n_markers, samples_alt_t(n_samples));
}

// ----------------------------------------------------------------------------
// Method population::simulate()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy> template <typename generator_t>
inline void population<n_ploidy>::simulate(generator_t && generator)
{
    simulate_founders_freqs(generator);
    simulate_founders_map(generator);

//    simulate_founders_alts(generator);
    fit_founders_alts();

    assign_samples_alts();
}

// ----------------------------------------------------------------------------
// Method population::simulate_founders_freqs()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy> template <typename generator_t>
inline void population<n_ploidy>::simulate_founders_freqs(generator_t && generator)
{
    SEQAN_ASSERT_EQ(founders_freqs.size(), n_founders);

    std::uniform_int_distribution<uint32_t> distribution(0u, n_founders - 1);
//    std::normal_distribution<uint32_t> distribution(0u, n_founders - 1);

    for (uint32_t i = 0; i < n_samples * n_ploidy - n_founders; i++)
        founders_freqs[distribution(generator)]++;

    std::cerr << "FOUNDERS FREQUENCIES: " << founders_freqs << std::endl;
    std::cerr << "=================================================================" << std::endl << std::endl;
}

// ----------------------------------------------------------------------------
// Method population::simulate_founders_map()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy> template <typename generator_t>
inline void population<n_ploidy>::simulate_founders_map(generator_t && generator)
{
    SEQAN_ASSERT_EQ(founders_map.num_elements(), n_samples * n_ploidy);
    SEQAN_ASSERT_EQ(std::accumulate(founders_freqs.begin(), founders_freqs.end(), 0u), n_samples * n_ploidy);

    auto founders_map_begin = founders_map.data();
    auto founders_map_it = founders_map.data();
    auto founders_map_end = founders_map.data() + founders_map.num_elements();
    for (auto founders_freqs_it = founders_freqs.begin(); founders_freqs_it != founders_freqs.end(); ++founders_freqs_it)
        founders_map_it = std::fill_n(founders_map_it, *founders_freqs_it, founders_freqs_it - founders_freqs.begin());
    SEQAN_ASSERT(founders_map_it == founders_map_end);

    std::shuffle(founders_map_begin, founders_map_end, generator);

    std::cerr << "FOUNDERS MAP: " << founders_map << std::endl;
    std::cerr << "=================================================================" << std::endl << std::endl;
}

// ----------------------------------------------------------------------------
// Method population::simulate_founders_alts()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy> template <typename generator_t>
inline void population<n_ploidy>::simulate_founders_alts(generator_t && generator)
{
    std::bernoulli_distribution distribution;

    // For each marker.
    std::for_each(founders_alts.begin(), founders_alts.end(), [&distribution, &generator](auto & founders_alt)
    {
        // Generate alleles for all founders.
        std::generate(founders_alt.begin(), founders_alt.end(), [&distribution, &generator]()
        {
            return distribution(generator);
        });
    });

    std::for_each(founders_alts.begin(), founders_alts.end(), [](auto & founders_alt)
    {
        std::cerr << "FOUNDERS ALLELE: " << founders_alt << std::endl;
    });
    std::cerr << "=================================================================" << std::endl << std::endl;
}


// ----------------------------------------------------------------------------
// Method population::fit_founders_alts()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
inline void population<n_ploidy>::fit_founders_alts()
{
    for (uint32_t marker_id = 0; marker_id < n_markers; ++marker_id)
        fit_founders_alts(marker_id);

    std::cerr << "=================================================================" << std::endl << std::endl;

    std::for_each(founders_alts.begin(), founders_alts.end(), [](auto & founders_alt)
    {
        std::cerr << "FOUNDERS ALLELE: " << founders_alt << std::endl;
    });
    std::cerr << "=================================================================" << std::endl << std::endl;
}

// ----------------------------------------------------------------------------
// Method population::fit_founders_alts()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
inline void population<n_ploidy>::fit_founders_alts(uint32_t marker_id)
{
    typedef typename dosages<lemon::Mip::Col, n_ploidy>::type mip_cols_t;

    lemon::Mip mip;
    mip_cols_t mip_z;                       // [dosage]
    mip_cols_t mip_dosages;                 // [dosage]
    std::vector<mip_cols_t> mip_errors;     // [sample][dosage]
    std::vector<mip_cols_t> mip_indicators; // [sample][dosage]
    std::vector<lemon::Mip::Col> mip_alts;  // [founder]

    auto const & marker_dosages = markers_in.dosages[marker_id];

    // Vector z to linearize l1-norm objective.
    std::generate(mip_z.begin(), mip_z.end(), [&mip]() { return mip.addCol(); });
    std::for_each(mip_z.begin(), mip_z.end(), [&mip](auto z)
    {
        mip.colType(z, lemon::Mip::REAL);
        mip.colLowerBound(z, 0);
    });

    // Dosages d.
    std::generate(mip_dosages.begin(), mip_dosages.end(), [&mip]() { return mip.addCol(); });
    std::for_each(mip_dosages.begin(), mip_dosages.end(), [&mip](auto d)
    {
        mip.colType(d, lemon::Mip::INTEGER);
        mip.colLowerBound(d, 0);
    });

    // Dosage absolute errors e_s,p.
    mip_errors.resize(n_samples);
    std::for_each(mip_errors.begin(), mip_errors.end(), [&mip](auto & mip_sample)
    {
        std::generate(mip_sample.begin(), mip_sample.end(), [&mip]() { return mip.addCol(); });
        std::for_each(mip_sample.begin(), mip_sample.end(), [&mip](auto e)
        {
            mip.colType(e, lemon::Mip::REAL);
            mip.colLowerBound(e, 0);
        });
    });

    // Dosage indicators i_s,p.
    mip_indicators.resize(n_samples);
    std::for_each(mip_indicators.begin(), mip_indicators.end(), [&mip](auto & mip_sample)
    {
        std::generate(mip_sample.begin(), mip_sample.end(), [&mip]() { return mip.addCol(); });
        std::for_each(mip_sample.begin(), mip_sample.end(), [&mip](auto i)
        {
            mip.colType(i, lemon::Mip::INTEGER);
            mip.colLowerBound(i, 0);
//            mip.colUpperBound(i, 1);
        });
    });

    // Founder alleles f.
    mip_alts.resize(n_founders);
    std::generate(mip_alts.begin(), mip_alts.end(), [&mip]() { return mip.addCol(); });
    std::for_each(mip_alts.begin(), mip_alts.end(), [&mip](auto a)
    {
        mip.colType(a, lemon::Mip::INTEGER);
        mip.colLowerBound(a, 0);
        mip.colUpperBound(a, 1);
    });

    // Constraint d - d_obs <= z.
    // Constraint d_obs - d <= z.
    for (uint8_t dosage = 0; dosage < n_ploidy + 1; dosage++)
    {
        mip.addRow(mip_dosages[dosage] - marker_dosages[dosage] <= mip_z[dosage]);
        mip.addRow(marker_dosages[dosage] - mip_dosages[dosage] <= mip_z[dosage]);
    }

    // Constraint d_p = Σ_s i_s,p.
    for (uint8_t dosage = 0; dosage < n_ploidy + 1; dosage++)
    {
        lemon::Mip::Expr i_sum;
        for (uint32_t sample_id = 0; sample_id < n_samples; sample_id++)
            i_sum += mip_indicators[sample_id][dosage];
        mip.addRow(mip_dosages[dosage] == i_sum);
    }

    // Constraint Σ_p i_s,p = 1.
    for (uint32_t sample_id = 0; sample_id < n_samples; sample_id++)
    {
        lemon::Mip::Expr i_sum;
        for (uint8_t dosage = 0; dosage < n_ploidy + 1; dosage++)
            i_sum += mip_indicators[sample_id][dosage];
        mip.addRow(i_sum == 1u);
    }

    // Constraint e_s,p >= Σ h(f) - p.
    // Constraint e_s,p >= p - Σ h(f).
    // Constraint i_s,p <= k * (1 - e_s,p).
    for (uint32_t sample_id = 0; sample_id < n_samples; sample_id++)
    {
        lemon::Mip::Expr sample_sum;
        for (uint8_t haplotype_id = 0; haplotype_id < n_ploidy; haplotype_id++)
            sample_sum += mip_alts[founders_map[sample_id][haplotype_id]];

        for (uint8_t dosage = 0; dosage < n_ploidy + 1; dosage++)
        {
            mip.addRow(mip_errors[sample_id][dosage] >= dosage - sample_sum);
            mip.addRow(mip_errors[sample_id][dosage] >= sample_sum - dosage);

            mip.addRow(mip_errors[sample_id][dosage] <= n_ploidy * (1u - mip_indicators[sample_id][dosage]));
        }
    }

    // Minimize Σ z.
    lemon::Mip::Expr o = std::accumulate(mip_z.begin(), mip_z.end(), lemon::Mip::Expr());
    mip.obj(o);
    mip.min();
//    mip.messageLevel(lemon::LpBase::MESSAGE_VERBOSE);
    mip.solve();

    SEQAN_ASSERT_EQ(mip.type(), lemon::MipSolver::OPTIMAL);

    typename dosages<uint32_t, n_ploidy>::type dosages_out;
    std::transform(mip_dosages.begin(), mip_dosages.end(), dosages_out.begin(), [&mip](auto d){ return mip.sol(d); });
    std::cerr << "DISTANCE: " << mip.solValue() << " = " << std::make_pair(marker_dosages, dosages_out) << std::endl;

//    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl << std::endl;
//
//    std::for_each(mip_errors.begin(), mip_errors.end(), [&mip](auto const & sample_errors)
//    {
//        dosages_t errors;
//        std::transform(sample_errors.begin(), sample_errors.end(), errors.begin(), [&mip](auto i){ return mip.sol(i); });
//        std::cerr << "ERRORS: " << errors << std::endl;
//    });
//
//    std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl << std::endl;
//
//    std::for_each(mip_indicators.begin(), mip_indicators.end(), [&mip](auto const & sample_indicators)
//    {
//        typename dosages<bool, n_ploidy>::type indicators;
//        std::transform(sample_indicators.begin(), sample_indicators.end(), indicators.begin(), [&mip](auto i){ return mip.sol(i); });
//        std::cerr << "INDICATORS: " << indicators << std::endl;
//    });
//
//    std::cerr << "-----------------------------------------------------------------" << std::endl << std::endl;

    // Update founder alleles.
    std::transform(mip_alts.begin(), mip_alts.end(), founders_alts[marker_id].begin(), [&mip](auto a){ return mip.sol(a); });
}

// ----------------------------------------------------------------------------
// Method population::assign_samples_alts()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
inline void population<n_ploidy>::assign_samples_alts()
{
    for (uint32_t marker_id = 0; marker_id < n_markers; ++marker_id)
        for (uint32_t sample_id = 0; sample_id < n_samples; ++sample_id)
            for (uint8_t haplotype_id = 0; haplotype_id < n_ploidy; ++haplotype_id)
                samples_alts[marker_id][sample_id][haplotype_id] = founders_alts[marker_id][founders_map[sample_id][haplotype_id]];

    std::for_each(samples_alts.begin(), samples_alts.end(), [](auto & samples_alt)
    {
        std::cerr << "SAMPLES ALLELE: " << samples_alt << std::endl;
    });
    std::cerr << "=================================================================" << std::endl << std::endl;

//    std::for_each(samples_alts.begin(), samples_alts.end(), [](auto & samples_alt)
//    {
//        std::cerr << "SAMPLES DOSAGES: " << alt_dosages<n_ploidy>(samples_alt) << std::endl;
//    });
//    std::cerr << "=================================================================" << std::endl << std::endl;
}

// ----------------------------------------------------------------------------
// Method population::write()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
inline void population<n_ploidy>::write(std::string const & /* vcf_filename_out */)
{
    // Open output file.
//    seqan::VcfFileOut vcf_out(seqan::toCString(options.vcf_filename_out));
    seqan::VcfFileOut vcf_out(std::cout, seqan::Vcf());

    // Fill contig names.
//    contigNames(context(vcf_out)) = contigNames(context(vcf_in));
    appendValue(contigNames(context(vcf_out)), "chr01");

    // Fill sample names.
    seqan::resize(sampleNames(context(vcf_out)), n_samples);
    for (uint32_t sample = 0; sample < n_samples; ++sample)
    {
        sampleNames(context(vcf_out))[sample] = "SAMPLE_";
        seqan::appendNumber(sampleNames(context(vcf_out))[sample], sample);
    }

    // Write VCF header.
    seqan::VcfHeader header_out;
    appendValue(header_out, seqan::VcfHeaderRecord("fileformat", "VCFv4.2"));
    appendValue(header_out, seqan::VcfHeaderRecord("ID", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    writeHeader(vcf_out, header_out);

    // Fill VCF record prototype.
    seqan::VcfRecord record_out;
    record_out.filter = "";
    record_out.info = "";
    record_out.format = "GT";
    seqan::resize(record_out.genotypeInfos, n_samples);

    // Write VCF records.
    for (uint32_t marker_id = 0; marker_id < n_markers; ++marker_id)
    {
        record_out.rID = std::get<0>(markers_in.positions[marker_id]);
        record_out.beginPos = std::get<1>(markers_in.positions[marker_id]);
        clear(record_out.id);
        appendNumber(record_out.id, marker_id);
        record_out.ref = front(markers_in.alleles[marker_id]);
        record_out.alt = back(markers_in.alleles[marker_id]);

        // Write VCF sample columns.
        for (uint32_t sample_id = 0; sample_id < n_samples; ++sample_id)
            write_genotype(record_out.genotypeInfos[sample_id], samples_alts[marker_id][sample_id]);

        writeRecord(vcf_out, record_out);
    }
}

// ============================================================================
// App
// ============================================================================

// ----------------------------------------------------------------------------
// Class app_options
// ----------------------------------------------------------------------------

struct app_options
{
    std::string vcf_filename_in;
    std::string vcf_filename_out;

    uint32_t seed;

    uint32_t n_ploidy;
    uint32_t n_founders;
    uint32_t n_samples;
    uint32_t n_markers;

    app_options():
        seed(0),
        n_ploidy(4),
        n_founders(1),
        n_samples(1),
        n_markers(1)
    {}

    void setup(seqan::ArgumentParser & parser) const
    {
        setAppName(parser, "SimBA-hap");
        setShortDescription(parser, "Haplotypes simulator");
        setCategory(parser, "Simulation");

        setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
        setDate(parser, SEQAN_DATE);

        addUsageLine(parser, "[\\fIOPTIONS\\fP]");

        addOption(parser, seqan::ArgParseOption("g", "seed", "Initial seed for pseudo-random number generation.",
                                                seqan::ArgParseOption::INTEGER));
        setMinValue(parser, "seed", "0");
        setDefaultValue(parser, "seed", seed);

        addOption(parser, seqan::ArgParseOption("v", "vcf", "Input VCF file.", seqan::ArgParseOption::INPUT_FILE));
        setValidValues(parser, "vcf", seqan::VcfFileIn::getFileExtensions());

        addOption(parser, seqan::ArgParseOption("p", "ploidy", "Organism ploidy.",
                                                seqan::ArgParseOption::INTEGER));
        setMinValue(parser, "ploidy", "2");
        setMaxValue(parser, "ploidy", "8");
        setDefaultValue(parser, "ploidy", n_ploidy);

        addOption(parser, seqan::ArgParseOption("f", "founders", "Number of founders to simulate.",
                                                seqan::ArgParseOption::INTEGER));
        setMinValue(parser, "founders", "1");
        setDefaultValue(parser, "founders", n_founders);

        addOption(parser, seqan::ArgParseOption("s", "samples", "Number of samples to simulate. Default: all samples in the input VCF file.",
                                                seqan::ArgParseOption::INTEGER));
        setMinValue(parser, "samples", "1");

        addOption(parser, seqan::ArgParseOption("m", "markers", "Number of markers to use. Default: all markers in the input VCF file.",
                                                seqan::ArgParseOption::INTEGER));
        setMinValue(parser, "markers", "1");
    }

    void parse(seqan::ArgumentParser const & parser)
    {
        getOptionValue(seed, parser, "seed");

        getOptionValue(vcf_filename_in, parser, "vcf");
        getOptionValue(n_ploidy, parser, "ploidy");
        getOptionValue(n_founders, parser, "founders");
        getOptionValue(n_samples, parser, "samples");
        getOptionValue(n_markers, parser, "markers");
    }
};

// ----------------------------------------------------------------------------
// Function run()
// ----------------------------------------------------------------------------

template <uint8_t n_ploidy>
void run(app_options const & options)
{
    std::mt19937 generator(options.seed);

    markers<n_ploidy> markers_in;
    population<n_ploidy> pop_out(markers_in);

    uint32_t n_markers = options.n_markers;

    if (options.vcf_filename_in.empty())
    {
        markers_in.simulate(options.n_samples, options.n_markers, generator);
    }
    else
    {
        markers_in.read(options.vcf_filename_in);
        n_markers = markers_in.dosages.size();
    }

    markers_in.normalize(options.n_samples);

    pop_out.resize(options.n_founders, options.n_samples, n_markers);
    pop_out.simulate(generator);

    pop_out.write(options.vcf_filename_out);
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    seqan::ArgumentParser parser;
    app_options options;

    options.setup(parser);

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.parse(parser);

    try
    {
        switch (options.n_ploidy)
        {
//            case 2:
//                run<2>(options);
//                break;
            case 4:
                run<4>(options);
                break;
//            case 6:
//                run<6>(options);
//                break;
//            case 8:
//                run<8>(options);
//                break;
            default:
                throw seqan::RuntimeError("Unsupported ploidy");
        }
    }
    catch (seqan::Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
