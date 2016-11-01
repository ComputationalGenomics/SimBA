// ==========================================================================
//                              SimBA-hap
// ==========================================================================
// Copyright (c) 2016 IBM Research

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>

//#include <boost/algorithm/string/join.hpp>

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    seqan::CharString vcf_filename_in;
    seqan::CharString vcf_filename_out;

    unsigned seed;

    unsigned n_ploidy;
    unsigned n_samples;
    unsigned n_crosses;

    bool input_phased;
    bool output_phased;

     Options():
        seed(0),
        n_ploidy(4),
        n_samples(std::numeric_limits<unsigned>::max()),
        n_crosses(1),
        input_phased(false),
        output_phased(false)
    {}
};

// ============================================================================
// Functions
// ============================================================================

//        for (auto const & c : contigNames(context(vcf_in)))
//            std::cerr << c << std::endl;
//        std::cerr << std::endl;
//        for (auto const & s : sampleNames(context(vcf_in)))
//            std::cerr << s << std::endl;
//        std::cerr << std::endl;
//
//        unsigned n = 10;
//        String<unsigned> zeros;
//        resize(zeros, n, 0, Exact());
//        StringSet<String<unsigned> > m;
//        for (unsigned i = 0; i < n; i++)
//            appendValue(m, zeros);
//        for (unsigned i = 0; i < n; i++)
//            for (unsigned j = i + 1; j < n; j++)
//                m[i][j] = i * n + j - (i + 1) * (i + 2) / 2;
//
//        for (unsigned i = 0; i < n; i++)
//        {
//            for (unsigned j = 0; j < n; j++)
//                std::cerr << m[i][j] << '\t';
//            std::cerr << std::endl;
//        }
//
//        for (unsigned x = 0; x < binom2(n); x++)
//        {
//            // Solve x = i * n + j - (i) * (i + 1) / 2 - 1 for i
//            double i_d = (2 * n - 1 - std::sqrt(4 * n * n - 4 * n - 8 * x - 7)) / 2;
//            unsigned i = (i_d == std::floor(i_d)) ? i_d - 1: i_d;
//            unsigned j = x - i * n + i * (i + 1) / 2 + i + 1;
//            std::cerr << x << " ==> " << std::make_pair(i,j) << std::endl;
//            std::cerr << binom2_idx_make_pair(x, n) << " <==> " << x << std::endl;
//        }

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
}

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

// ----------------------------------------------------------------------------
// Function run()
// ----------------------------------------------------------------------------

void run(Options const & options)
{
    typedef std::pair<unsigned, unsigned> samples_comb_t;
    typedef std::pair<unsigned, unsigned> haplos_comb_t;
    typedef std::tuple<samples_comb_t, haplos_comb_t, haplos_comb_t> f1_comb_t;
    typedef std::vector<char> genotype_t;

    // Open input file.
    seqan::VcfFileIn vcf_in(seqan::toCString(options.vcf_filename_in));

    // Read header.
    seqan::VcfHeader header_in;
    seqan::readHeader(header_in, vcf_in);
    auto sample_names = sampleNames(context(vcf_in));
    auto n_samples = std::min(options.n_samples, static_cast<unsigned>(seqan::length(sample_names)));

    // Generate a random F1 population from input samples.
    std::vector<f1_comb_t> f1s;
    f1s.reserve(options.n_crosses);
    std::mt19937 generator(options.seed);
    rand_f1(n_samples, options.n_ploidy, options.n_crosses, f1s, generator);

    // Open output file.
//    seqan::VcfFileOut vcf_out(seqan::toCString(options.vcf_filename_out));
    seqan::VcfFileOut vcf_out(std::cout, seqan::Vcf());
    contigNames(context(vcf_out)) = contigNames(context(vcf_in));

    seqan::VcfHeader header_out;
    appendValue(header_out, seqan::VcfHeaderRecord("fileformat", "VCFv4.2"));
    appendValue(header_out, seqan::VcfHeaderRecord("ID", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
//    appendValue(header_out, seqan::VcfHeaderRecord("phasing", "partial"));

    genotype_t genotype;
    std::vector<genotype_t> sample_genotypes;
    std::vector<genotype_t> f1_genotypes;
    sample_genotypes.reserve(n_samples);
    f1_genotypes.reserve(options.n_crosses);

    // Write F1 cross names.
//    for (auto const & f1 : f1s)
//    {
//        auto parents_f1 = std::get<0>(f1);
//        auto haplos_p1 = std::get<1>(f1);
//        auto haplos_p2 = std::get<2>(f1);

//        std::string cross_name = sample_names[std::get<0>(parents_f1)] + "H" + std::get<0>(haplos_p1);
//
//            sample_names[std::get<0>(parents_f1)][std::get<1>(haplos_p1)] + "|" +
//            sample_names[std::get<1>(parents_f1)][std::get<0>(haplos_p2)] + "|" +
//            sample_names[std::get<1>(parents_f1)][std::get<1>(haplos_p2)])

//    appendValue(sampleNames(context(vcf_out)), "NA00001");

//        appendValue(sampleNames(context(vcf_out)),
//            sample_names[std::get<0>(parents_f1)] [std::get<0>(haplos_p1)] + "|" +
//            sample_names[std::get<0>(parents_f1)][std::get<1>(haplos_p1)] + "|" +
//            sample_names[std::get<1>(parents_f1)][std::get<0>(haplos_p2)] + "|" +
//            sample_names[std::get<1>(parents_f1)][std::get<1>(haplos_p2)])
//    }

//    writeHeader(vcf_out, header_out);

    for (auto const & f1 : f1s)
        std::cerr << "INFO: " << f1 << std::endl;

    std::mt19937 phaser(options.seed);
    std::mt19937 shuffler(options.seed);

    seqan::VcfRecord record_in;
    seqan::VcfRecord record_out;

//    std::string genotype_out;
//    const char * const genotype_delim = "|";

    // Read the file recordwise.
    while (!atEnd(vcf_in))
    {
        sample_genotypes.clear();
        f1_genotypes.clear();

        readRecord(record_in, vcf_in);

        // Read all samples.
        for (auto const & genotype_info : record_in.genotypeInfos)
        {
            genotype.clear();

            // Read genotype string.
            auto genotype_it = directionIterator(genotype_info, seqan::Input());
            seqan::readUntil(genotype, genotype_it, seqan::EqualsChar<':'>());

            // Convert genotype string to vector by removing slashes.
            genotype.erase(std::remove(genotype.begin(), genotype.end(), '/'), genotype.end());

            // Unknown genotypes are taken as hom ref.
            if (genotype.size() == 1 & genotype[0] == '.')
                genotype = {'0', '0', '0', '0'};

            if (genotype.size() != options.n_ploidy)
                throw seqan::ParseError("Wrong ploidy.");

            sample_genotypes.push_back(genotype);
        }

        if (!options.input_phased)
        {
            // Convert genotypes to canonical form.
//                std::sort(genotype.begin(), genotype.end(), std::greater<char>());

            // Simulate random phases.
            for (auto & genotype : sample_genotypes)
                std::shuffle(genotype.begin(), genotype.end(), phaser);
        }

        // Compose F1s.
        for (auto const & f1 : f1s)
        {
            auto parents_f1 = std::get<0>(f1);
            auto haplos_p1 = std::get<1>(f1);
            auto haplos_p2 = std::get<2>(f1);

            f1_genotypes.push_back({
                sample_genotypes[std::get<0>(parents_f1)][std::get<0>(haplos_p1)],
                sample_genotypes[std::get<0>(parents_f1)][std::get<1>(haplos_p1)],
                sample_genotypes[std::get<1>(parents_f1)][std::get<0>(haplos_p2)],
                sample_genotypes[std::get<1>(parents_f1)][std::get<1>(haplos_p2)]
            });
        }

        // Shuffle F1 phases.
        if (!options.output_phased)
        {
            for (auto & genotype : f1_genotypes)
                std::shuffle(genotype.begin(), genotype.end(), shuffler);
        }

        // Write all F1 crosses.
        record_out.rID = record_in.rID;
        record_out.beginPos = record_in.beginPos;
        record_out.id = record_in.id;
        record_out.ref = record_in.ref;
        record_out.alt = record_in.alt;
        record_out.filter = record_in.filter;
        record_out.info = "";
        record_out.format = "GT";
//        for (auto const & genotype : f1_genotypes)
//        {
//            auto genotype_info = boost::algorithm::join(genotype, "|");
//            appendValue(record_out.genotypeInfos, genotype_info);
//        }

//        writeRecord(vcf_out, record_out);

        for (auto const & genotype : f1_genotypes)
            std::copy(genotype.begin(), genotype.end(), std::ostream_iterator<char>(std::cout, "\t"));
        std::cout << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(seqan::ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "vcf_sim");
    setShortDescription(parser, "VCF Simulator");
    setCategory(parser, "Simulation");

    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIVCF FILE IN\\fP> <\\fIVCF FILE OUT\\fP>");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE));
    setValidValues(parser, 0, seqan::VcfFileIn::getFileExtensions());

//    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_FILE));
//    setValidValues(parser, 1, seqan::VcfFileOut::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("g", "seed", "Initial seed for pseudo-random number generation.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "seed", "0");
    setDefaultValue(parser, "seed", options.seed);

    addOption(parser, seqan::ArgParseOption("s", "samples", "Number of input samples to consider. Default: all samples in the input VCF file.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "samples", "2");
    setMaxValue(parser, "samples", "1000");

    addOption(parser, seqan::ArgParseOption("c", "crosses", "Number of F1 crosses to simulate.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "crosses", "1");
    setMaxValue(parser, "crosses", "1000");
    setDefaultValue(parser, "crosses", options.n_crosses);

    addOption(parser, seqan::ArgParseOption("ip", "input-phased", "Input is already phased."));
    addOption(parser, seqan::ArgParseOption("op", "output-phased", "Output is kept phased."));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, seqan::ArgumentParser & parser, int argc, char const ** argv)
{
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Parse vcf input filename.
    getArgumentValue(options.vcf_filename_in, parser, 0);
//    getArgumentValue(options.vcf_filename_out, parser, 1);

    getOptionValue(options.seed, parser, "seed");
    getOptionValue(options.n_samples, parser, "samples");
    getOptionValue(options.n_crosses, parser, "crosses");

    getOptionValue(options.input_phased, parser, "input-phased");
    getOptionValue(options.output_phased, parser, "output-phased");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    seqan::ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    try
    {
        run(options);
    }
    catch (seqan::Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

