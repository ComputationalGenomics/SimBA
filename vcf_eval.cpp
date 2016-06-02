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

