//
//  main.cpp
//  tk_textile
//
//  Created by Joel Davis on 2/18/16.
//  Copyright © 2016 Joel Davis. All rights reserved.
//

#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define TK_OBJFILE_IMPLEMENTATION
#include "tk_objfile.h"

#include "optionparser.h"
#include "tk_textiler.h"

enum OptionIndex {
    Opt_UNKNOWN,

    Opt_Help,
    
    Opt_SourceImageFilename,
    Opt_SourceMeshFilename,
    Opt_DestTextureFilename,
    
    Opt_EdgeSize,
};

static option::ArgStatus ArgRequired(const option::Option& option, bool msg)
{
    if (option.arg != 0)
        return option::ARG_OK;
    
    if (msg) printf("Option '%s' requires an argument\n", option.name );
    
    return option::ARG_ILLEGAL;
}

static option::ArgStatus ArgNumeric(const option::Option& option, bool msg)
{
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
        return option::ARG_OK;
    
    if (msg) printf("Option '%s' requires a numeric argument\n", option.name );
    return option::ARG_ILLEGAL;
}


const option::Descriptor usage[] = {
    { Opt_UNKNOWN, 0, "", "",  option::Arg::None, "Usage: blah blah blah..." },
    
    { Opt_Help, 0, "", "help", option::Arg::None, "  --help     \tPrint this help message." },
    { Opt_SourceImageFilename, 0, "i", "image", ArgRequired, "  --image, -i \tSource Image filename." },
    { Opt_SourceMeshFilename, 0, "m", "mesh", ArgRequired, "  --mesh, -m \tSource Mesh .OBJ filename." },
    { Opt_DestTextureFilename, 0, "ot", "outtex", ArgRequired, "  --outtex, -ot \tOutput texture filename." },
    { Opt_EdgeSize, 0, "e", "edgesize", ArgNumeric, "  --edgesize, -e \tOutput edge size." },
    
    {0,0,0,0,0,0}
};

int main(int argc, const char * argv[])
{
    tapnik::TextureTiler *textiler = new tapnik::TextureTiler();
    
    argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
    option::Stats stats(usage, argc, argv);
    
    option::Option *options = new option::Option[ stats.options_max ];
    option::Option *buffer = new option::Option[ stats.buffer_max];
 
    option::Parser parse(usage, argc, argv, options, buffer);
    
    if (parse.error())
        return 1;

    if (options[Opt_Help] || argc == 0)
    {
        int columns = getenv("COLUMNS")? atoi(getenv("COLUMNS")) : 80;
        option::printUsage(fwrite, stdout, usage, columns);
        return 0;
    }
    
    const char *outTexFilename = "outtexture.png";
    
    for (int i = 0; i < parse.optionsCount(); ++i)
    {
        option::Option& opt = buffer[i];
        switch (opt.index())
        {
            case Opt_SourceImageFilename:
                printf("Source image is %s\n", opt.arg);
                textiler->sourceImage_ = tapnik::Image::load( opt.arg );                
                break;
            case Opt_SourceMeshFilename:
                printf("Source mesh is %s\n", opt.arg);
                textiler->mesh_ = tapnik::Mesh::load( opt.arg );
                break;
            case Opt_EdgeSize:
                printf("Edge size is '%s'\n", opt.arg );
                textiler->edgeSize_ = atoi( opt.arg );
                break;
            case Opt_DestTextureFilename:
                outTexFilename = opt.arg;
                break;
        }
    }
    
    delete [] options;
    delete [] buffer;
    
    if (!textiler->sourceImage_) {
        printf("Must specify a source image: --image, --i <image>\n");
    } else {
        // Actually do stuff
        textiler->outTexFilename_ = strdup(outTexFilename);
        textiler->doStuff( outTexFilename );
    }
    
    delete textiler;
}
