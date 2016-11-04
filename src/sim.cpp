// ==========================================================================
//                              SimBA-hap
// ==========================================================================
// Copyright (c) 2016 IBM Research

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>

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

//template <typename T1, typename T2, typename T3>
//std::ostream& operator<<(std::ostream & s, std::tuple<T1, T2, T3> const & t)
//{
//    s << '{' << std::get<0>(t) << "; " << std::get<1>(t) << "; " << std::get<2>(t) << '}';
//    return s;
//}

template <typename T>
std::ostream& operator<<(std::ostream & s, std::vector<T> const & v)
{
    s << "[ ";
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(s, " "));
    s << "]";
    return s;
}
}

// ============================================================================
// Functions on genotypes
// ============================================================================

// ----------------------------------------------------------------------------
// Function read_alleles()
// ----------------------------------------------------------------------------
// Read genotype string.

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

    auto genotype_it = directionIterator(genotype_info, seqan::Input());
    seqan::readUntil(genotype, genotype_it, seqan::EqualsChar<':'>());

    // Convert genotype string to vector by removing slashes.
    genotype.erase(std::remove(genotype.begin(), genotype.end(), '/'), genotype.end());
    genotype.erase(std::remove(genotype.begin(), genotype.end(), '|'), genotype.end());
}

// ----------------------------------------------------------------------------
// Function ref_genotype()
// ----------------------------------------------------------------------------

//template <typename genotype_t>
//inline void ref_genotype(genotype_t & genotype)
//{
//    std::fill(genotype.begin(), genotype.end(), 0u);
//}

// ----------------------------------------------------------------------------
// Function get_ploidy()
// ----------------------------------------------------------------------------

template <typename genotype_t>
inline uint32_t get_ploidy(genotype_t const & genotype)
{
    return genotype.size();
}

// ----------------------------------------------------------------------------
// Function ref_dosage()
// ----------------------------------------------------------------------------

template <typename genotype_t>
inline size_t ref_dosage(genotype_t const & genotype)
{
    return std::count(genotype.begin(), genotype.end(), '0');
}

// ----------------------------------------------------------------------------
// Function is_unknown()
// ----------------------------------------------------------------------------

template <typename genotype_t>
inline bool is_unknown(genotype_t const & genotype)
{
    return genotype.size() == 1u && genotype.front() == '.';
}

// ----------------------------------------------------------------------------
// Function is_homozygous()
// ----------------------------------------------------------------------------

//template <typename genotype_t>
//inline bool is_homozygous(genotype_t const & genotype)
//{
//    return std::adjacent_find(genotype.begin(), genotype.end(), std::not_equal_to<char>()) == genotype.end();
//}

// ============================================================================
// Markers
// ============================================================================

struct markers
{
    typedef std::pair<uint32_t, uint32_t> position_t;
//    typedef std::array<seqan::IupacString, 2> alleles_t;
    typedef seqan::StringSet<seqan::IupacString, seqan::Owner<seqan::ConcatDirect<>>> alleles_t;
    typedef std::vector<uint32_t> dosages_t;

    std::vector<position_t> positions; // positions[m] = (chr, pos)
    std::vector<alleles_t> alleles;    // alleles[m] = (ref, alt)
    std::vector<dosages_t> dosages;    // dosages[m] = [f(0),f(1),...,f(p)]
};

// ----------------------------------------------------------------------------
// Method read_markers()
// ----------------------------------------------------------------------------

inline void read_markers(markers & markers_in, uint32_t n_ploidy, std::string const & vcf_filename_in)
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
            std::cerr << "VARIANT @ " << position << " POLYALLELIC" << std::endl;
            continue;
        }

        dosages.resize(n_ploidy);
        std::fill(dosages.begin(), dosages.end(), 0u);

        // Read all sample genotypes.
        for (auto const & genotype_info : vcf_record.genotypeInfos)
        {
            read_genotype(genotype, genotype_info);

            // Skip unknown genotypes.
            if (is_unknown(genotype))
            {
                std::cerr << "GENOTYPE @ " << position << " UNKNOWN" << std::endl;
                continue;
            }

            // Stop on wrong ploidy.
            if (get_ploidy(genotype) != n_ploidy)
                throw seqan::IOError("Input ploidy does not match VCF genotypes");

            dosages[ref_dosage(genotype)]++;
        }

        markers_in.positions.push_back(position);
        markers_in.alleles.push_back(alleles);
        markers_in.dosages.push_back(dosages);

        std::cerr << "DOSAGES @ " << position << " # " << dosages << std::endl;
    }

    std::cerr << "=================================================================" << std::endl << std::endl;
}

// ============================================================================
// Founders
// ============================================================================

typedef std::vector<uint32_t> founders;

// ----------------------------------------------------------------------------
// Method simulate_founders_freqs()
// ----------------------------------------------------------------------------

template <typename founders_freq_t, typename generator_t>
inline void simulate_founders_freqs(founders_freq_t & founders_freqs, uint32_t n_founders, uint32_t n_samples, uint32_t n_ploidy, generator_t && generator)
{
    SEQAN_ASSERT_GT(n_samples, 0u);
    SEQAN_ASSERT_GEQ(n_samples, n_founders);

    std::uniform_int_distribution<uint32_t> distribution(1u, n_founders - 1);
//    std::normal_distribution<uint32_t> distribution(1u, n_founders - 1);

    founders_freqs.resize(n_founders, 1u);

    for (uint32_t i = 0; i < n_samples * n_ploidy - n_founders; i++)
        founders_freqs[distribution(generator)]++;

    std::cerr << "FOUNDERS: " << founders_freqs << std::endl;
    std::cerr << "=================================================================" << std::endl << std::endl;
}

// ----------------------------------------------------------------------------
// Method simulate_founders_assignment()
// ----------------------------------------------------------------------------

template <typename haplotypes_t, typename founders_freq_t, typename generator_t>
inline void simulate_founders_assignment(haplotypes_t & haplotypes, founders_freq_t const & founders_freqs, uint32_t n_samples, uint32_t n_ploidy, generator_t && generator)
{
    SEQAN_ASSERT_EQ(std::accumulate(founders_freqs.begin(), founders_freqs.end(), 0u), n_samples * n_ploidy);

    haplotypes.resize(n_ploidy * n_samples);

    auto haplotypes_it = haplotypes.begin();
    for (auto founders_freqs_it = founders_freqs.begin(); founders_freqs_it != founders_freqs.end(); ++founders_freqs_it)
        haplotypes_it = std::fill_n(haplotypes_it, *founders_freqs_it, founders_freqs_it - founders_freqs.begin());
    SEQAN_ASSERT(haplotypes_it == haplotypes.end());

    std::shuffle(haplotypes.begin(), haplotypes.end(), generator);

    std::cerr << "HAPLOTYPES: " << haplotypes << std::endl;
    std::cerr << "=================================================================" << std::endl << std::endl;
}

// ============================================================================
// Population
// ============================================================================

struct population
{
    typedef std::vector<bool> marker_t;

    std::vector<marker_t> founders;                 // founders[m][f] == ALT
//    std::vector<std::vector<uint32_t>> haplotypes;  // haplotypes[n][p] = i
    std::vector<uint32_t> haplotypes;
};

// ============================================================================
// App options
// ============================================================================

// ----------------------------------------------------------------------------
// Class simba_hap_options
// ----------------------------------------------------------------------------

class simba_hap_options
{
public:
    std::string vcf_filename_in;

    uint32_t seed;

    uint32_t n_ploidy;
    uint32_t n_founders;
    uint32_t n_samples;
    uint32_t n_markers;

    void parse(seqan::ArgumentParser const & parser)
    {
        getOptionValue(seed, parser, "seed");

        getOptionValue(vcf_filename_in, parser, "vcf");
        getOptionValue(n_ploidy, parser, "ploidy");
        getOptionValue(n_founders, parser, "founders");
        getOptionValue(n_samples, parser, "samples");
        getOptionValue(n_markers, parser, "markers");
    }

    simba_hap_options():
        seed(0),
        n_ploidy(4),
        n_founders(1),
        n_samples(1),
        n_markers(std::numeric_limits<uint32_t>::max())
    {}
};

// ----------------------------------------------------------------------------
// Function setup_argument_parser()
// ----------------------------------------------------------------------------

void setup_argument_parser(seqan::ArgumentParser & parser, simba_hap_options const & options)
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
    setDefaultValue(parser, "seed", options.seed);

    addOption(parser, seqan::ArgParseOption("v", "vcf", "Input VCF file.", seqan::ArgParseOption::INPUT_FILE));
    setValidValues(parser, "vcf", seqan::VcfFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("p", "ploidy", "Organism ploidy.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "ploidy", "2");
    setMaxValue(parser, "ploidy", "8");
    setDefaultValue(parser, "ploidy", options.n_ploidy);

    addOption(parser, seqan::ArgParseOption("f", "founders", "Number of founders to simulate.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "founders", "1");
    setDefaultValue(parser, "founders", options.n_founders);

    addOption(parser, seqan::ArgParseOption("s", "samples", "Number of samples to simulate. Default: all samples in the input VCF file.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "samples", "1");

    addOption(parser, seqan::ArgParseOption("m", "markers", "Number of markers to use. Default: all markers in the input VCF file.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "markers", "1");
}

// ============================================================================
// App
// ============================================================================

// ----------------------------------------------------------------------------
// Class simba_hap
// ----------------------------------------------------------------------------

struct simba_hap
{
    // Input.
    founders founders_in;
    markers markers_in;

    // Output.
    population pop_out;

    void run(simba_hap_options const & options);
};

// ----------------------------------------------------------------------------
// Method simba_hap::run()
// ----------------------------------------------------------------------------

void simba_hap::run(simba_hap_options const & options)
{
    std::mt19937 generator(options.seed);

    simulate_founders_freqs(founders_in, options.n_founders, options.n_samples, options.n_ploidy, generator);
    simulate_founders_assignment(pop_out.haplotypes, founders_in, options.n_samples, options.n_ploidy, generator);

    if (!options.vcf_filename_in.empty())
        read_markers(markers_in, options.n_ploidy, options.vcf_filename_in);
//    else
//        simulate_markers(markers_in, options.n_ploidy, options.n_samples, options.n_markers);
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    seqan::ArgumentParser parser;
    simba_hap app;
    simba_hap_options options;

    setup_argument_parser(parser, options);
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.parse(parser);

    try
    {
        app.run(options);
    }
    catch (seqan::Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
