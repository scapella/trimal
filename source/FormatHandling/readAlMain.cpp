/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

#include "FormatHandling/FormatManager.h"
#include "Alignment/Alignment.h"
#include "defines.h"
#include "values.h"

#include <iostream>
#include <cstring>
#include <iomanip>

int parseArguments(int argc, char *argv[],
        FormatHandling::FormatManager* machine,
        std::vector<std::string>* inFiles,
        std::vector<std::string>* outFormats,
        std::string* outPattern)
{
    if (argc == 1) return 1;
    for(int i = 1; i < argc; i++ )
    {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
            return -1;
        // Check if current argument is the '-in' argument.
        if (!strcmp(argv[i], "-in"))
        {
            // Check if this is the last argument.
            if (i >= argc -1)
            {
                std::cerr << "ERROR: At least one file should be passed after the '-in' argument\n";
                return 1;
            }
            
            // Check if the next argument is a parameter or file argument.
            else if (argv[i + 1][0] == '-')
            {
                std::cerr << "ERROR: At least one file should be passed after the '-in' argument and you passed argument " << argv[i + 1] << "\n";
                return 1;
            }
            
            // Add every file that is not a parameter
            else while(++i != argc)
            {
                if (argv[i][0] == '-')
                {
                    i--;
                    break;
                }
                else
                {
                    inFiles->emplace_back(argv[i]);
                }
            }
        }
        
        // Check if current argument is the '-out' argument.
        else if (!strcmp(argv[i], "-out"))
        {
            // Check if this is the last argument.
            if (i >= argc -1)
            {
                std::cerr << "A file pattern should be passed after the '-out' argument\n";
                return 1;
            }
            else 
            {
                if (argv[i + 1][0] == '-')
                {
                    std::cerr << "A file pattern should be passed after the '-out' argument and you passed argument " << argv[i + 1] << "\n";
                    return 1;
                }
                *outPattern = argv[++i];
                continue;
            }
        }
        
        // Check if current argument is the '-formats' argument.
        else if (!strcmp(argv[i], "-formats"))
        {
            if (i >= argc -1)
            {
                std::cerr << "A format should be passed after the '-formats' argument\n";
                return 1;
            }
            else while(++i != argc)
            {
                if (argv[i][0] == '-')
                {
                    i--;
                    break;
                }
                else
                {
                    outFormats->emplace_back(argv[i]);
                }
            }
        }
        
        else if (!strcmp(argv[i], "-reverse"))
        {
            machine->reverse = true;
        }
        else if (!strcmp(argv[i], "-keepHeaders"))
        {
            machine->keepHeader = true;
        }
        
        //Compatibility with legacy options:

        else if (!strcmp(argv[i], "-html"))
            outFormats->emplace_back("html");
        
        else if (!strcmp(argv[i], "-nbrf"))
            outFormats->emplace_back("nbrf");
        
        else if (!strcmp(argv[i], "-mega"))
            outFormats->emplace_back("mega");
        
        else if (!strcmp(argv[i], "-nexus"))
            outFormats->emplace_back("nexus");
        
        else if (!strcmp(argv[i], "-clustal"))
            outFormats->emplace_back("clustal");
        
        else if (!strcmp(argv[i], "-fasta") || !strcmp(argv[i], "-onlyseqs"))
            outFormats->emplace_back("fasta");
        
        else if (!strcmp(argv[i], "-fasta_m10"))
        {
            outFormats->emplace_back("fasta");
        }
        
        else if (!strcmp(argv[i], "-phylip"))
            outFormats->emplace_back("phylip40");
        
        else if (!strcmp(argv[i], "-phylip_m10"))
        {
            outFormats->emplace_back("phylip40_m10");
        }
        
        else if (!strcmp(argv[i], "-phylip_paml"))
            outFormats->emplace_back("phylippaml");
        
        else if (!strcmp(argv[i], "-phylip_paml_m10"))
        {
            outFormats->emplace_back("phylippaml_m10");
        }
        
        else if (!strcmp(argv[i], "-phylip3.2"))
            outFormats->emplace_back("phylip32");
        
        else if (!strcmp(argv[i], "-phylip3.2_m10"))
        {
            outFormats->emplace_back("phylip32_m10");
        }
        else if (!strcmp(argv[i], "-format"))
        {
            machine->format = true;
        }
        else if (!strcmp(argv[i], "-type"))
        {
            machine->type = true;
        }
        else if (!strcmp(argv[i], "-info"))
        {
            machine->info = true;
        }

        // If a command is not recognized, give an error.
        else
        {
            std::cerr << argv[i] << " not recognized or repeated.\n";
            return 1;
        }
    }

#if debug
    std::cout << "Input Files:\n";
    for (std::string ifile : *inFiles)
    {
        std::cout << "-> Input file: " << ifile << "\n";
    }
    
    std::cout << "Out Formats:\n";
    for (std::string oformat : *outFormats)
    {
        std::cout << "-> Output format: " << oformat << "\n";
    }
    
    std::cout << "Under the pattern\n-> " << *outPattern << "\n";
#endif
    
    return 0;
}

int checkArguments(FormatHandling::FormatManager* machine, std::vector<std::string>* inFiles, std::vector<std::string>* outFormats, std::string* outPattern)
{
    int returnValue = 0;
    if (inFiles->size() == 0)
    {
        std::cerr        << "ERROR: At least one input file must be provided\n";
        returnValue = 1;
    }
    if (*outPattern == "")
    {
        if (inFiles->size() == 1 && outFormats->size() == 1 && !(machine->format || machine->info || machine->type))
            machine->hasOutputFile = false;
        else if (outFormats->size() != 0)
        {
            std::cerr    << "ERROR: Terminal output option not compatible with information printing (-info | -format | -type)\n" 
                    << "Provide an output format or disable information printing.\n";
            returnValue = 1;
        }
    }
    else if (outFormats->size() == 0)
    {
        std::cerr        << "ERROR: At least one output format must be provided\n";
        returnValue = 1;
    }
    return returnValue;
}

void menu()
{

    std::cout << "\n"
    << "readAl v" << VERSION << ".rev" << REVISION << " build[" << BUILD
    << "]. " << AUTHORS << "\n\n"

    << "readAl webpage: http://trimal.cgenomics.org\n\n"

    << "This program is free software: you can redistribute it and/or modify "
    << "\n"
    << "it under the terms of the GNU General Public License as published by "
    << "\n"
    << "the Free Software Foundation, the last available version.\n"
    << "\n"

    << "BASIC USAGE\n\n"
    << "\treadalMS -in <inputfiles> -out <pattern> -format [formats] [options].\n\n"

    << "\t-h                    Show this information.\n"
//     << "\t--version            Show readAl version.\n\n"

    << "\t-in <inputfiles>      Input files in several formats. Separated by spaces.\n"
    << "\t                      Available formats are: " << FormatHandling::FormatManager().getInputFormatsAvailable() << "\n" 
    << "\t-out <pattern>        Output file name pattern (default STDOUT).\n"
    << "\t                      It will replace optional the tags [in]        -> Original filename without extension.\n"
    << "\t                                                        [format]    -> Output's format name\n"
    << "\t                                                        [extension] -> Output's extension\n"
    << "\n"
    
    << "\t-formats             Formats you want the output to be converted to.\n"
    << "\t                     Available formats are: " << FormatHandling::FormatManager().getOutputFormatsAvailable() << "\n" 
    << "\t                     Being the HTML format not a format itself, but a colored report of the alignment files.\n\n"
    << "\t-format              Print information about input file format "
    << "and if sequences are aligned or not.\n"

    << "\t-type                Print information about biological "
    << "sequences datatype (e.g. nucleotides:dna, nucleotides:rna, aminoacids, etc)"
    << "\n"

    << "\t-info                Print information about sequences number, "
    << "average sequence length, max & min sequence length"
    << "\n" 
    
    << "\t-reverse             Output the reverse of sequences in "
    << "input file.\n\n"

    << "\t-keepHeaders         Keeps the headers of the original format if it had any\n\n"

    
    << "LEGACY OPTIONS\nTake in mind that this arguments may be discontinued any time."<< "\n\n"
    
    << "\t-onlyseqs            Generate output with only residues from "
    << "input file\n\n"
    
    << "\t-html                Output residues colored according their "
    << "physicochemical properties. HTML file.\n\n"


    << "\t-nbrf                Output file in NBRF/PIR format\n"
    << "\t-mega                Output file in MEGA format\n"

    << "\t-nexus               Output file in NEXUS format\n"
    << "\t-clustal             Output file in CLUSTAL format\n"
    << "\n"

    << "\t-fasta               Output file in FASTA format\n"
    << "\t-fasta_m10           Output file in FASTA format. Sequences "
    << "name up to 10 characters.\n\n"

    << "\t-phylip              Output file in PHYLIP/PHYLIP4 format"
    << "\n"
    << "\t-phylip_m10          Output file in PHYLIP/PHYLIP4 format. "
    << "Sequences name up to 10 characters.\n"
    << "\t-phylip_paml         Output file in PHYLIP format compatible "
    << "with PAML\n"
    << "\t-phylip_paml_m10     Output file in PHYLIP format compatible "
    << "with PAML. Sequences name up to 10 characters.\n"
    << "\t-phylip3.2           Output file in PHYLIP3.2 format\n"
    << "\t-phylip3.2_m10       Output file in PHYLIP3.2 format. Sequences"
    << " name up to 10 characters.\n\n"
    << "If you specify any m10 format, this will result in all formats having the sequences names shortened as this has the same effect as '-shortNames' argument\n\n"
    
    
    << "EXAMPLES OF USE\n\n"
    
    << "\treadalMS -in ./dataset/AA1.fas -out ./dataset/[in].output.[extension] -formats clustal\n"
    << "\t  -> Will produce ./dataset/AA1.output.clw\n\n"
    
    << "\treadalMS -in ./dataset/example1.clw -out ./dataset/[in].[format].[extension] -formats fasta phylip32 phylip40\n"
    << "\t  -> Will produce ./dataset/example1.FASTA.fasta ./dataset/example1.PHYLIP32.phy ./dataset/example1.PHYLIP40.phy\n\n"
    
    << "\treadalMS -in ./dataset/example1.clw -out ./dataset/[in]/[format].[extension] -formats fasta phylip32 phylip40\n"
    << "\t  -> Will produce ./dataset/example1/FASTA.fasta ./dataset/example1/PHYLIP32.phy ./dataset/example1/PHYLIP40.phy\n"
    << "\t     ONLY if ./dataset/example1/ already exists.\n\n" 
    
    << "\treadalMS -in ./dataset/AA1.fas ./dataset/AA2.fas -out ./dataset/[in].output.[extension] -formats clustal pir\n"
    << "\t  -> Will produce ./dataset/AA1.output.clw ./dataset/AA2.output.clw ./dataset/AA1.output.pir ./dataset/AA2.output.pir\n\n"
    
    << "\treadalMS -in ./dataset/AA1.fas -format -type -info\n"
    << "\t  -> Will produce terminal output giving information about AA1.fas alignment file\n\n"
    
    << "\treadalMS -in ./dataset/AA1.fas ./dataset/AA2.fas -out ./dataset/[in].output.[extension] -formats html\n"
    << "\t  -> Will produce ./dataset/AA1.output.html ./dataset/AA2.output.html\n"
    << "\t     Those files are not indeed reformats of the original alignments, but an HTML colored report of the alignment file.\n"
    
    << "\n";
}


int main(int argc, char *argv[])
{
    FormatHandling::FormatManager MachineState = FormatHandling::FormatManager();
    
    std::vector<std::string> outFormats = std::vector<std::string>();
    std::vector<std::string> inFiles  = std::vector<std::string>();
    std::string outPattern;
    
    int result = parseArguments(argc, argv, &MachineState, &inFiles, &outFormats, &outPattern);
    if (result == 1)
    {
        menu();
        return 0;
    }
    else if (result != 0) return result;

    result = checkArguments(&MachineState, &inFiles, &outFormats, &outPattern);
    if (result != 0) return result;

    if(MachineState.format || MachineState.info || MachineState.type) {
        for (const std::string &str : inFiles)
        {
            Alignment* alignment = MachineState.loadAlignment(str);
            if (alignment != nullptr)
            {

                std::cout << "## Alignment File:\t" << str << "\n";
                
                if (MachineState.format)
                    /* Inform about if sequences are aligned or not */
                    std::cout    << "## Input file format\t" << MachineState.getFileFormatName(str) << "\n"
                            << "## Input file aligned\t" << (alignment->isAligned ? "YES":"NO")
                            << "\n";

                if(MachineState.type) {
                    /* Inform about biological datatype */
                    if (alignment->getAlignmentType() == SequenceTypes::DNA)
                        std::cout << "## Input file datatype\tnucleotides:dna\n";
                    else if (alignment->getAlignmentType() == (SequenceTypes::DNA | SequenceTypes::DEG))
                        std::cout << "## Input file datatype\tnucleotides:dna_degenerate_codes\n";
                    else if (alignment->getAlignmentType() == SequenceTypes::RNA)
                        std::cout << "## Input file datatype\tnucleotides:rna\n";
                    else if (alignment->getAlignmentType() == (SequenceTypes::RNA | SequenceTypes::DEG))
                        std::cout << "## Input file datatype\tnucleotides:rna_degenerate_codes\n";
                    else if (alignment->getAlignmentType() == SequenceTypes::AA)
                        std::cout << "## Input file datatype\tamino-acids\n";
                    else if (alignment->getAlignmentType() == (SequenceTypes::AA | SequenceTypes::DEG))
                        std::cout << "## Input file datatype\tamino-acids_degenerate_codes\n";
                    else
                        std::cout << "## Input file datatype\tunknown\n";
                }

                if(MachineState.info)
                    alignment->printAlignmentInfo(std::cout);
                
                std::cout << "\n";
                
                if (!outFormats.empty())
                    MachineState.saveAlignment(outPattern, outFormats, *alignment);
                
            }
            delete alignment;
        }
    }
    else if (!outFormats.empty() || MachineState.reverse)
        MachineState.loadAndSaveMultipleAlignments(inFiles, outPattern, outFormats);
    else
        std::cerr << "ERROR: An option has to be chosen\n";
}
