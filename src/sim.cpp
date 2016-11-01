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

template <typename T1, typename T2, typename T3>
std::ostream& operator<<(std::ostream & s, std::tuple<T1, T2, T3> const & t)
{
    s << '{' << std::get<0>(t) << "; " << std::get<1>(t) << "; " << std::get<2>(t) << '}';
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
}

// ============================================================================
// Functions on genotypes
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
// Function is_homozygous()
// ----------------------------------------------------------------------------

template <typename genotype_t>
inline bool is_homozygous(genotype_t const & genotype)
{
    return std::adjacent_find(genotype.begin(), genotype.end(), std::not_equal_to<char>()) == genotype.end();
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

template <typename genotype_t>
inline void ref_genotype(genotype_t & genotype)
{
    std::fill(genotype.begin(), genotype.end(), 0u);
}

// ----------------------------------------------------------------------------
// Function read_alleles()
// ----------------------------------------------------------------------------
// Read genotype string.

template <typename alleles_t>
inline void read_alleles(alleles_t & alleles, seqan::VcfRecord const & record)
{
    front(alleles) = record.ref;

    clear(back(alleles));
    auto alleles_it = directionIterator(record.alt, seqan::Input());
    seqan::readUntil(back(alleles), alleles_it, seqan::EqualsChar<','>());
}

// ============================================================================
// Markers
// ============================================================================

struct markers
{
    typedef std::pair<uint32_t, uint32_t> position_t;
    typedef std::array<seqan::IupacString, 2> alleles_t;
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

            // Skip unknown variants.
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

// ----------------------------------------------------------------------------
// Method simulate_founders_freq()
// ----------------------------------------------------------------------------

template <typename founders_freq_t>
inline void simulate_founders_freq(founders_freq_t & founders_freqs_in, uint32_t n_founders)
{

}

// ============================================================================
// Functions for F1 simulation
// ============================================================================

// ----------------------------------------------------------------------------
// Function binom2()
// ----------------------------------------------------------------------------

template <typename T>
inline T binom2(T n)
{
    return n * (n - 1) / 2;
}

// ----------------------------------------------------------------------------
// Function binom2_idx_make_pair()
// ----------------------------------------------------------------------------

template <typename T, typename T2>
inline std::pair<T,T> binom2_idx_make_pair(T x, T2 n)
{
    // Solve x = i * n + j - i * (i + 1) / 2 - 1 for i.
    double i_d = (2 * n - 1 - std::sqrt(4 * n * n - 4 * n - 8 * x - 7)) / 2;

    T i = (i_d == std::floor(i_d)) ? i_d - 1: i_d;
    T j = x - i * n + i * (i + 1) / 2 + i + 1;

    return std::make_pair(i, j);
}

// ----------------------------------------------------------------------------
// Function rand_f1()
// ----------------------------------------------------------------------------

template <typename T1, typename T2, typename T3, typename F1, typename URNG>
void rand_f1(T1 n_samples, T2 n_ploidy, T3 n_crosses, F1 & f1s, URNG && g)
{
    // Number of distinct F1 crosses.
    auto n_parent_comb = binom2(n_samples);
    auto n_haplo_comb = binom2(n_ploidy);
    auto n_meiosis_comb = n_haplo_comb * n_haplo_comb;
    auto n_f1_comb = n_parent_comb * n_meiosis_comb;

    // Generate all random F1 crosses.
    std::vector<uint64_t> crosses_ids(n_f1_comb);
    std::iota(crosses_ids.begin(), crosses_ids.end(), 0);
    std::shuffle(crosses_ids.begin(), crosses_ids.end(), g);

    // Select first n crosses.
    for (auto crosses_it = crosses_ids.begin(); crosses_it != crosses_ids.begin() + n_crosses; ++crosses_it)
    {
        auto cross_id = *crosses_it;
        auto id_parents_f1 = cross_id / n_meiosis_comb;
        auto id_haplos_f1 = cross_id % n_meiosis_comb;
        auto id_haplos_p1 = id_haplos_f1 % n_haplo_comb;
        auto id_haplos_p2 = id_haplos_f1 / n_haplo_comb;

        auto parents_f1 = binom2_idx_make_pair(id_parents_f1, n_samples);
        auto haplos_p1 = binom2_idx_make_pair(id_haplos_p1, n_ploidy);
        auto haplos_p2 = binom2_idx_make_pair(id_haplos_p2, n_ploidy);

        SEQAN_ASSERT_LT(parents_f1.first, n_samples);
        SEQAN_ASSERT_LT(parents_f1.second, n_samples);
        SEQAN_ASSERT_LT(haplos_p1.first, n_ploidy);
        SEQAN_ASSERT_LT(haplos_p1.second, n_ploidy);
        SEQAN_ASSERT_LT(haplos_p2.first, n_ploidy);
        SEQAN_ASSERT_LT(haplos_p2.second, n_ploidy);

        f1s.push_back(std::make_tuple(parents_f1, haplos_p1, haplos_p2));
    }
}

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
    uint32_t n_samples;
    uint32_t n_variants;

    void parse(seqan::ArgumentParser const & parser)
    {
        getOptionValue(seed, parser, "seed");

        getOptionValue(vcf_filename_in, parser, "vcf");
        getOptionValue(n_ploidy, parser, "ploidy");
        getOptionValue(n_samples, parser, "samples");
        getOptionValue(n_variants, parser, "variants");
    }

     simba_hap_options():
        seed(0),
        n_ploidy(4),
        n_samples(std::numeric_limits<uint32_t>::max()),
        n_variants(std::numeric_limits<uint32_t>::max())
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

    addOption(parser, seqan::ArgParseOption("f", "vcf", "Input VCF file.", seqan::ArgParseOption::INPUT_FILE));
    setValidValues(parser, "vcf", seqan::VcfFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("p", "ploidy", "Organism ploidy.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "ploidy", "2");
    setMaxValue(parser, "ploidy", "8");

    addOption(parser, seqan::ArgParseOption("s", "samples", "Number of input samples to use. Default: all samples in the input VCF file.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "samples", "1");
    setMaxValue(parser, "samples", "1000");

    addOption(parser, seqan::ArgParseOption("v", "variants", "Number of variants to use. Default: all variants in the input VCF file.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "variants", "1");
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
    std::vector<uint32_t> founders_freqs_in;
    markers markers_in;

    // Output.
    std::vector<std::vector<bool>> founder_haplotypes_out; // founder_haplotypes_out[m] = [ h_i[m] == alt ].
    std::vector<std::vector<uint32_t>> sample_haplotypes_out; // sample_haplotypes_out[n][p] = i

    void run(simba_hap_options const & options);
};

// ----------------------------------------------------------------------------
// Method simba_hap::run()
// ----------------------------------------------------------------------------

void simba_hap::run(simba_hap_options const & options)
{
//    simulate_founders_freqs(founders_freqs_in, options.n_founders);

    if (!options.vcf_filename_in.empty())
        read_markers(markers_in, options.n_ploidy, options.vcf_filename_in);
//    else
//        simulate_markers(markers_in, options.n_ploidy, options.n_samples, options.n_variants);
}


//void simba_hap::run(simba_hap_options const & options)
//{
//    typedef std::pair<unsigned, unsigned> samples_comb_t;
//    typedef std::pair<unsigned, unsigned> haplos_comb_t;
//    typedef std::tuple<samples_comb_t, haplos_comb_t, haplos_comb_t> f1_comb_t;
//    typedef std::vector<char> genotype_t;
//
//    // Open input file.
//    seqan::VcfFileIn vcf_in(seqan::toCString(options.vcf_filename_in));
//
//    // Read header.
//    seqan::VcfHeader header_in;
//    seqan::readHeader(header_in, vcf_in);
//    auto sample_names = sampleNames(context(vcf_in));
//    auto n_samples = std::min(options.n_samples, static_cast<unsigned>(seqan::length(sample_names)));
//
//    // Generate a random F1 population from input samples.
//    std::vector<f1_comb_t> f1s;
//    f1s.reserve(options.n_crosses);
//    std::mt19937 generator(options.seed);
//    rand_f1(n_samples, options.n_ploidy, options.n_crosses, f1s, generator);
//
//    genotype_t genotype;
//    std::vector<genotype_t> sample_genotypes;
//    std::vector<genotype_t> f1_genotypes;
//    sample_genotypes.reserve(n_samples);
//    f1_genotypes.reserve(options.n_crosses);
//
//    for (auto const & f1 : f1s)
//        std::cerr << "INFO: " << f1 << std::endl;
//
//    std::mt19937 phaser(options.seed);
//    std::mt19937 shuffler(options.seed);
//
//    seqan::VcfRecord record_in;
//    seqan::VcfRecord record_out;
//
//    // Read the file recordwise.
//    while (!atEnd(vcf_in))
//    {
//        sample_genotypes.clear();
//        f1_genotypes.clear();
//
//        readRecord(record_in, vcf_in);
//
//        // Read all samples.
//        for (auto const & genotype_info : record_in.genotypeInfos)
//        {
//            genotype.clear();
//
//            // Read genotype string.
//            auto genotype_it = directionIterator(genotype_info, seqan::Input());
//            seqan::readUntil(genotype, genotype_it, seqan::EqualsChar<':'>());
//
//            // Convert genotype string to vector by removing slashes.
//            genotype.erase(std::remove(genotype.begin(), genotype.end(), '/'), genotype.end());
//
//            // Unknown genotypes are taken as hom ref.
//            if (is_unknown(genotype))
//                genotype = {'0', '0', '0', '0'};
//
//            if (genotype.size() != options.n_ploidy)
//                throw seqan::ParseError("Wrong ploidy.");
//
//            sample_genotypes.push_back(genotype);
//        }
//
//        if (!options.input_phased)
//        {
//            // Convert genotypes to canonical form.
//            std::sort(genotype.begin(), genotype.end(), std::greater<char>());
//
//            // Simulate random phases.
//            for (auto & genotype : sample_genotypes)
//                std::shuffle(genotype.begin(), genotype.end(), phaser);
//        }
//
//        // Compose F1s.
//        for (auto const & f1 : f1s)
//        {
//            auto parents_f1 = std::get<0>(f1);
//            auto haplos_p1 = std::get<1>(f1);
//            auto haplos_p2 = std::get<2>(f1);
//
//            f1_genotypes.push_back({
//                sample_genotypes[std::get<0>(parents_f1)][std::get<0>(haplos_p1)],
//                sample_genotypes[std::get<0>(parents_f1)][std::get<1>(haplos_p1)],
//                sample_genotypes[std::get<1>(parents_f1)][std::get<0>(haplos_p2)],
//                sample_genotypes[std::get<1>(parents_f1)][std::get<1>(haplos_p2)]
//            });
//        }
//
//        // Shuffle F1 phases.
//        if (!options.output_phased)
//        {
//            for (auto & genotype : f1_genotypes)
//                std::shuffle(genotype.begin(), genotype.end(), shuffler);
//        }
//
//        // Output F1 genotypes.
//        for (auto const & genotype : f1_genotypes)
//            std::copy(genotype.begin(), genotype.end(), std::ostream_iterator<char>(std::cout, "\t"));
//        std::cout << std::endl;
//    }
//}

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
