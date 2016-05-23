// ==========================================================================
//                              VCF Simulator
// ==========================================================================
// Copyright (c) 2016, Enrico Siragusa, IBM Research
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Enrico Siragusa or IBM Research nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR IBM RESEARCH BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <esiragu@us.ibm.com>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    seqan::CharString vcf_filename_in;

    unsigned n_ploidy;
    unsigned n_crosses;

    bool input_phased;
    bool output_phased;

     Options():
        n_ploidy(4),
        n_crosses(5),
        input_phased(false),
        output_phased(false)
    {}
};

// ============================================================================
// Functions
// ============================================================================

//        for (auto const & c : contigNames(context(vcfIn)))
//            std::cerr << c << std::endl;
//        std::cerr << std::endl;
//        for (auto const & s : sampleNames(context(vcfIn)))
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

template <typename T1, typename T2, typename T3, typename F1>
void rand_f1(T1 n_samples, T2 n_ploidy, T3 n_crosses, F1 & f1s)
{
    // Number of distinct F1 crosses.
    auto n_parent_comb = binom2(n_samples);
    auto n_haplo_comb = binom2(n_ploidy);
    auto n_meiosis_comb = n_haplo_comb * n_haplo_comb;
    auto n_f1_comb = n_parent_comb * n_meiosis_comb;

    // Generate all random F1 crosses.
    std::vector<uint64_t> crosses_ids(n_f1_comb);
    std::iota(crosses_ids.begin(), crosses_ids.end(), 0);
    std::mt19937 generator;
    std::shuffle(crosses_ids.begin(), crosses_ids.end(), generator);

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
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(seqan::ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "vcf_sim");
    setShortDescription(parser, "VCF Simulator");
    setCategory(parser, "Simulation");

    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIVCF FILE IN\\fP>");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE));
    setValidValues(parser, 0, seqan::VcfFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("c", "crosses", "Number of F1 crosses to simulate.",
                                            seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "crosses", "0");
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

    getOptionValue(options.n_crosses, parser, "crosses");

    getOptionValue(options.input_phased, parser, "input-phased");
    getOptionValue(options.output_phased, parser, "output-phased");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function mainWithOptions()
// ----------------------------------------------------------------------------

int mainWithOptions(Options const & options)
{
    typedef std::pair<unsigned, unsigned> samples_comb_t;
    typedef std::pair<unsigned, unsigned> haplos_comb_t;
    typedef std::tuple<samples_comb_t, haplos_comb_t, haplos_comb_t> f1_comb_t;
    typedef std::vector<char> genotype_t;

    try
    {
        // Open input file.
        seqan::VcfFileIn vcfIn(seqan::toCString(options.vcf_filename_in));

        // Copy over header.
        seqan::VcfHeader header;
        seqan::readHeader(header, vcfIn);
        auto n_samples = seqan::length(sampleNames(context(vcfIn)));

        genotype_t genotype;
        std::vector<genotype_t> sample_genotypes;
        std::vector<genotype_t> f1_genotypes;
        sample_genotypes.reserve(n_samples);
        f1_genotypes.reserve(options.n_crosses);

        // Generate a random F1 population from input samples.
        std::vector<f1_comb_t> f1s;
        f1s.reserve(options.n_crosses);
        rand_f1(n_samples, options.n_ploidy, options.n_crosses, f1s);

        for (auto const & f1 : f1s)
            std::cerr << "INFO: " << f1 << std::endl;

//        std::cerr << "INFO: " << 0 << " distinct parent haplotypes" << std::endl;

        std::mt19937 phaser;

        // Read the file recordwise.
        seqan::VcfRecord record;
        while (!atEnd(vcfIn))
        {
            sample_genotypes.clear();
            f1_genotypes.clear();

            readRecord(record, vcfIn);

            // Read all samples.
            for (auto const & genotype_info : record.genotypeInfos)
            {
                genotype.clear();

                // Read genotype string.
                auto genotype_it = directionIterator(genotype_info, seqan::Input());
                seqan::readUntil(genotype, genotype_it, seqan::EqualsChar<':'>());

                // Convert genotype string to vector by removing /.
                genotype.erase(std::remove(genotype.begin(), genotype.end(), '/'), genotype.end());

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

                // DEBUG
//                auto genotype_p1 = sample_genotypes[std::get<0>(parents_f1)];
//                auto genotype_p2 = sample_genotypes[std::get<1>(parents_f1)];
//                std::cerr << f1 << "\t==>\t";
//                std::copy(genotype_p1.begin(), genotype_p1.end(), std::ostream_iterator<char>(std::cerr, " "));
//                std::cerr << "\tX\t";
//                std::copy(genotype_p2.begin(), genotype_p2.end(), std::ostream_iterator<char>(std::cerr, " "));
//                std::cerr << "\t=\t";
//                std::copy(f1_genotypes.back().begin(), f1_genotypes.back().end(), std::ostream_iterator<char>(std::cerr, " "));
//                std::cerr << std::endl;
            }

            // Shuffle F1 phases.
            if (!options.output_phased)
            {
                for (auto & genotype : f1_genotypes)
                    std::shuffle(genotype.begin(), genotype.end(), phaser);
            }

            // Write all F1 crosses.
            for (auto const & genotype : f1_genotypes)
            {
                std::copy(genotype.begin(), genotype.end(), std::ostream_iterator<char>(std::cout, "\t"));
//                std::cout << '\t';
            }
            std::cout << std::endl;
        }
    }
    catch (seqan::Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
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

    if (res == seqan::ArgumentParser::PARSE_OK)
        return mainWithOptions(options);
    else
        return res;
}

