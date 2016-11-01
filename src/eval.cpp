// ==========================================================================
//                              SimBA-hap
// ==========================================================================
// Copyright (c) 2016 IBM Research

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>

#include <boost/algorithm/string/join.hpp>

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    seqan::CharString vcf_filename_in;

    Options() {}
};

// ============================================================================
// Functions
// ============================================================================

template <typename P1, typename P2>
inline unsigned distance(P1 const & a, P2 & const b)
{

}

// ----------------------------------------------------------------------------
// Function run()
// ----------------------------------------------------------------------------

void run(Options const & /* options */)
{
    typedef std::vector<char> haplotype_t;
    typedef std::array<haplotype_t, 4> phase_t;

    phase_t solution = { {
        {'1', '1', '1', '0', '0', '0', '1'},
        {'1', '0', '1', '0', '0', '1', '0'},
        {'0', '0', '0', '1', '1', '0', '1'}
    } };

    phase_t candidate = { {
        {'1', '1', '1', '0', '0', '0', '0'},
        {'1', '0', '1', '1', '1', '0', '1'},
        {'0', '0', '0', '0', '0', '1', '1'}
    } };


}

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(seqan::ArgumentParser & parser, Options const & /* options */)
{
    setAppName(parser, "phase_eval");
    setShortDescription(parser, "VCF Simulator");

    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIVCF FILE IN\\fP>");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE));
    setValidValues(parser, 0, seqan::VcfFileIn::getFileExtensions());
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

