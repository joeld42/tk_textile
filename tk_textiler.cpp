//
//  tk_textile.cpp
//  tk_textile
//
//  Created by Joel Davis on 2/27/16.
//  Copyright © 2016 Joel Davis. All rights reserved.
//

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "tk_textiler.h"
#include "stb_image.h"
#include "stb_image_write.h"

#define TK_SQRT3_OVER_2 (0.86602540378443864676372)
#define TK_SGN(x) ((x<0)?-1:((x>0)?1:0))
#define TK_ABS(x) (((x)<0)?-(x):(x))

using namespace tapnik;

// ==================================================
# pragma mark - Misc
// ==================================================

float randUniform()
{
    return (float)rand() / (float)RAND_MAX;
}

float randUniform( float minVal, float maxVal )
{
    return minVal + (randUniform() * (maxVal-minVal));
}

void *readEntireFile( const char *filename, size_t *out_filesz )
{
    FILE *fp = fopen( filename, "r" );
    if (!fp) return NULL;
    
    // Get file size
    fseek( fp, 0L, SEEK_END );
    size_t filesz = ftell(fp);
    fseek( fp, 0L, SEEK_SET );
    
    void *fileData = malloc( filesz );
    if (fileData)
    {
        size_t result = fread( fileData, filesz, 1, fp );
        
        // result is # of chunks read, we're asking for 1, fread
        // won't return partial reads, so it's all or nothing.
        if (!result)
        {
            free( fileData);
            fileData = NULL;
        }
        else
        {
            // read suceeded, set out filesize
            *out_filesz = filesz;
        }
    }
    
    return fileData;
}

void loadObjErrorMessage( size_t lineNum, const char *message, void *userData )
{
    Mesh *mesh = (Mesh*)userData;
    printf("Error loading OBJ file '%s' on line %zu: %s\n",
           mesh?mesh->filename_:"",
           lineNum, message );
}

// ==================================================
# pragma mark - Image
// ==================================================

Image::Image( uint32_t width, uint32_t height ) :
    width_(width), height_(height)
{
    imgdata_ = (uint32_t*)malloc( sizeof(uint32_t)*4*width_*height_ );
}

Image::~Image()
{
    if (imgdata_) {
        free( imgdata_ );
    }
    
    if (filename_) {
        free( filename_ );
    }
}

void Image::clear( uint32_t color )
{
    // huh this is a neat function. i wonder if it's portable??
    memset_pattern4( imgdata_, &color, sizeof(uint32_t)*4*width_*height_ );
}

inline void Image::drawPixel( int32_t x, int32_t y, uint32_t color  )
{
    if ((x>=0) && (y>=0) && (x < width_) && (y < height_)) {
        imgdata_[ (y * width_) + x ] = color;
    }
}

inline uint32_t Image::getPixel( int32_t x, int32_t y )
{
    if ((x>=0) && (y>=0) && (x < width_) && (y < height_)) {
        return imgdata_[ (y * width_) + x ];
    } else {
        return 0xffff00ff; // noticable out of bounds color
    }
}

// shamelessly stolen off the internet
void Image::drawLine( int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color )
{
    int32_t dx=x2-x1;      /* the horizontal distance of the line */
    int32_t dy=y2-y1;      /* the vertical distance of the line */
    int32_t dxabs=TK_ABS(dx);
    int32_t dyabs=TK_ABS(dy);
    int32_t sdx=TK_SGN(dx);
    int32_t sdy=TK_SGN(dy);
    int32_t x=dyabs>>1;
    int32_t y=dxabs>>1;
    int32_t px=x1;
    int32_t py=y1;
    
    drawPixel(px,py,color );
    
    if (dxabs>=dyabs) { /* the line is more horizontal than vertical */
        for(int i=0;i<dxabs;i++) {
            y+=dyabs;
            if (y>=dxabs) {
                y-=dxabs;
                py+=sdy;
            }
            px+=sdx;
            drawPixel(px,py,color );
        }
    } else { /* the line is more vertical than horizontal */
        for(int i=0;i<dyabs;i++) {
            x+=dxabs;
            if (x>=dyabs)
            {
                x-=dyabs;
                px+=sdx;
            }
            py+=sdy;
            drawPixel(px,py,color );
        }
    }
}

void Image::drawFatLine( int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color )
{
    drawLine( x1, y1, x2, y2, color );
    drawLine( x1-1, y1, x2-1, y2, color );
    drawLine( x1+1, y1, x2+1, y2, color );
    drawLine( x1, y1-1, x2, y2-1, color );
    drawLine( x1, y1+1, x2, y2+1, color );
}
Image *Image::load( const char *filename )
{
    Image *img = new Image();
    img->filename_ = strdup(filename);
    
    int comp;
    img->imgdata_ = (uint32_t*)stbi_load( filename, (int*)&img->width_, (int*)&img->height_, &comp, 4 );
    
    if (!img->imgdata_) {
        printf("Error loading Image %s\n", filename );
        delete img;
        img = nullptr;
    } else {
        printf("Loaded Image %s (%dx%d)\n", filename, img->width_, img->height_ );
    }
    
    return img;
}

void Image::drawImage( int32_t x, int32_t y, tapnik::Image *img )
{
    for (int32_t j = 0; j < img->height_; j++) {
        for (int32_t i = 0; i < img->width_; i++) {
            // TODO: alpha blend. 
            drawPixel( x+i, y+j, img->getPixel(i, j));
        }
    }
        
}

void Image::save()
{
    if (stbi_write_png(filename_, width_, height_, 4, imgdata_, width_*4 )) {
        printf("Wrote output image: %s\n", filename_ );
    } else {
        printf("Error writing output image: %s\n", filename_);
    }
}

// ==================================================
# pragma mark - Mesh
// ==================================================


void meshProcessTriangle( TK_TriangleVert a, TK_TriangleVert b, TK_TriangleVert c, void *userData )
{
    Mesh *mesh = (Mesh*)userData;
    Triangle *t = mesh->meshTris_ + mesh->numMeshTris_++;
    t->A_ = a;
    t->B_ = b;
    t->C_ = c;
}

Mesh *Mesh::load(const char *filename)
{
    Mesh *result = new Mesh();
    
    // Callbacks for API
    TK_ObjDelegate objDelegate = {};
    objDelegate.userData = (void*)result;
    objDelegate.error = loadObjErrorMessage;
    
    size_t objFileSize = 0;
    void *objFileData = readEntireFile( filename, &objFileSize );
    if (!objFileData) {
        printf("Could not open .OBJ file '%s'\n", filename );
    }
    
    // Prepass to determine memory reqs
    TK_ParseObj( objFileData, objFileSize, &objDelegate );
    printf("Scratch Mem: %zu\n", objDelegate.scratchMemSize );
    objDelegate.scratchMem = malloc( objDelegate.scratchMemSize );
    
    result->numMeshTris_ = 0;
    result->meshTris_ = new Triangle[ objDelegate.numTriangles ];
    
    // Parse again with memory
    objDelegate.triangle = meshProcessTriangle;
    
    TK_ParseObj( objFileData, objFileSize, &objDelegate );

    printf("Num Triangles %zu\n",result->numMeshTris_);
    
    return result;
}

Mesh::~Mesh()
{
    if (filename_) {
        free(filename_);
    }
}

static inline bool vertMatch(TK_TriangleVert a, TK_TriangleVert b, float threshold = 0.001 )
{
    return ( (fabs(a.pos[0] - b.pos[0]) < threshold) &&
             (fabs(a.pos[1] - b.pos[1]) < threshold) &&
             (fabs(a.pos[2] - b.pos[2]) < threshold) );
}

static inline bool edgeMatch( TK_TriangleVert a1, TK_TriangleVert b1,
                TK_TriangleVert a2, TK_TriangleVert b2 )
{
    return ( (vertMatch(a1, a2) && vertMatch(b1, b2)) ||
             (vertMatch(a1, b2) && vertMatch(b1, a2)) );
}

static inline void checkEdge(TK_TriangleVert a1, TK_TriangleVert b1,
                             TK_TriangleVert a2, TK_TriangleVert b2,
                             Triangle *triA, Triangle *triB,
                             Triangle **targA, Triangle **targB ) {
    if ((!*targA)&&(edgeMatch(a1, b1, a2, b2))) {
        *targA = triB;
        *targB = triA;
    }
}

void Mesh::buildAdjacency()
{
   printf("build adjacency...\n");
    for (Triangle *tri = meshTris_; (tri - meshTris_) < numMeshTris_; tri++) {
        for (Triangle *other = meshTris_; (other - meshTris_) < numMeshTris_; other++) {
            
            if (tri==other) continue;
            
            // If we found all of our neighbors, we're done
            if ((tri->nbAB_) && (tri->nbBC_) && (tri->nbCA_)) break;
            
            // Look for neighbor across AB
            checkEdge( tri->A_, tri->B_, other->A_, other->B_, tri, other, &(tri->nbAB_), &(other->nbAB_));
            checkEdge( tri->A_, tri->B_, other->B_, other->C_, tri, other, &(tri->nbAB_), &(other->nbBC_));
            checkEdge( tri->A_, tri->B_, other->C_, other->A_, tri, other, &(tri->nbAB_), &(other->nbCA_));

            // Look for neighbor across BC
            checkEdge( tri->B_, tri->C_, other->A_, other->B_, tri, other, &(tri->nbBC_), &(other->nbAB_));
            checkEdge( tri->B_, tri->C_, other->B_, other->C_, tri, other, &(tri->nbBC_), &(other->nbBC_));
            checkEdge( tri->B_, tri->C_, other->C_, other->A_, tri, other, &(tri->nbBC_), &(other->nbCA_));

            // Look for neighbor across CA
            checkEdge( tri->C_, tri->A_, other->A_, other->B_, tri, other, &(tri->nbCA_), &(other->nbAB_));
            checkEdge( tri->C_, tri->A_, other->B_, other->C_, tri, other, &(tri->nbCA_), &(other->nbBC_));
            checkEdge( tri->C_, tri->A_, other->C_, other->A_, tri, other, &(tri->nbCA_), &(other->nbCA_));
        }
    }
    
    printf("build adjacency done...\n");
}

static void assignBackEdge( Triangle *tri, Triangle *nbr, EdgeInfo *edge )
{
    if (nbr->nbAB_ == tri) {
        nbr->ab_ = edge;
    } else if (nbr->nbBC_ == tri) {
        nbr->bc_ = edge;
    } else if (nbr->nbCA_ == tri) {
        nbr->ca_ = edge;
    } else {
        printf("WARN: couldn't find back edge for nbr.\n");
    }
}

void Mesh::assignEdges()
{
    printf("Assign edges\n");
    
    EdgeInfo *edges[3];
    edges[0] = new EdgeInfo( 1, 0xffff0000 );
    edges[1] = new EdgeInfo( 2, 0xff00ff00 );
    edges[2] = new EdgeInfo( 3, 0xff0000ff );
    
    for (Triangle *tri = meshTris_; (tri - meshTris_) < numMeshTris_; tri++) {
        // for now, just pick at random
        if (!tri->ab_) {
            tri->ab_ = edges[ rand() % 3 ];
            assignBackEdge( tri, tri->nbAB_, tri->ab_ );
        }
        
        if (!tri->bc_) {
            tri->bc_ = edges[ rand() % 3 ];
            assignBackEdge( tri, tri->nbBC_, tri->bc_ );
        }
        
        if (!tri->ca_) {
            tri->ca_ = edges[ rand() % 3 ];
            assignBackEdge( tri, tri->nbCA_, tri->ca_ );
        }
    }
    
    printf("Assign edges done...\n");
}

// ==================================================
# pragma mark - Tile
// ==================================================
Tile::Tile( int32_t edgeSize, int32_t margin ) {
    
    img_ = new tapnik::Image( edgeSize + (margin*2), (uint32_t)((float)edgeSize * TK_SQRT3_OVER_2) + (margin*2) );

    tileA_[0] = margin + (img_->width_ / 2);
    tileA_[1] = margin;
    
    tileB_[0] = (img_->width_) - margin;
    tileB_[1] = (img_->height_-1) - margin;
    
    tileC_[0] = margin;
    tileC_[1] = (img_->height_-1) - margin;
    
    xform_ = GLKMatrix4Identity;
}

void Tile::debugDrawAnnotations()
{
    img_->drawFatLine(tileA_[0], tileA_[1], tileB_[0], tileB_[1], edge_[0]->debugColor_ );
    img_->drawFatLine(tileB_[0], tileB_[1], tileC_[0], tileC_[1], edge_[1]->debugColor_ );
    img_->drawFatLine(tileC_[0], tileC_[1], tileA_[0], tileA_[1], edge_[2]->debugColor_ );
}

void Tile::paintFromSource(Image *srcImage)
{
    // TMP: DEBUG assign a random transform
    xform_ = GLKMatrix4MakeTranslation( randUniform( 0, srcImage->width_ - img_->width_),
                                        randUniform( 0, srcImage->height_ - img_->height_),
                                        0.0 );
    
    for (int j=0; j < img_->height_; j++) {
        for (int i=0; i < img_->width_; i++) {
            
            // TODO: test point for including in triangle
            
            GLKVector4 samplePos = GLKVector4Make( (float)i, (float)j, 0.0, 1.0 );
            samplePos = GLKMatrix4MultiplyVector4( xform_, samplePos);
            
            // TODO: (maybe) fractional lookup and interpolate
            uint32_t sampleVal = srcImage->getPixel( (int32_t)samplePos.x, (int32_t)samplePos.y );
            img_->drawPixel(i, j, sampleVal );
        }
    }
}

// ==================================================
# pragma mark - TextureTiler
// ==================================================
TextureTiler::~TextureTiler()
{
    delete sourceImage_;
}


// Does stuff.
void TextureTiler::doStuff( const char *outTexFilename )
{
    // build mesh
    mesh_->buildAdjacency();
    mesh_->assignEdges();
    
    // make tiles
    gatherTiles();
    
    debugDumpTiles();
    
    assembleTiles();

//    // TMP: copy a chunk of src image into the out image
//    for (int j=0; j < outTexture_->height_; j++) {
//        memcpy( outTexture_->imgdata_ + (outTexture_->width_)*j,
//                sourceImage_->imgdata_ + (sourceImage_->width_)*j,
//                outTexture_->width_ * 4 );
//    }
    
    finish();
}

Tile *TextureTiler::findOrCreateTile( Triangle *tri )
{
    // Note: could be clever here and support rotated tiles, e.g. use tile ABC for
    // BCA but i'm not sure the number of tiles is an issue, and this will give more
    // variety in the image anyways so not doing it yet...
    Tile *result = nullptr;
    for (size_t tileNdx = 0; tileNdx < numTiles_; tileNdx++) {
        Tile *tile = tiles_[tileNdx];
        if ( (tri->ab_ == tile->edge_[0]) &&
             (tri->bc_ == tile->edge_[1]) &&
             (tri->ca_ == tile->edge_[2]) ) {
            result = tile;
            break;
        }
    }
    
    // Didn't find one, create a new tile
    if (!result) {
        // Do we need space for more tiles?
        if (numTiles_ == tilesCapacity_) {
            size_t growSize = (tilesCapacity_ * 3) / 2;
            if (growSize < 10) growSize = 10;
            
            tilesCapacity_ += growSize;
            if (!tiles_) {
                tiles_ = (Tile**)malloc(tilesCapacity_ * sizeof(Tile*));
            } else {
                tiles_ = (Tile**)realloc( tiles_, tilesCapacity_ * sizeof(Tile*));
            }
            
            printf("Grow tiles, size %zu capacity %zu\n", numTiles_, tilesCapacity_ );
        }
        
        result = new Tile( edgeSize_, 4 );
        result->edge_[0] = tri->ab_;
        result->edge_[1] = tri->bc_;
        result->edge_[2] = tri->ca_;
        tiles_[numTiles_++] = result;
    }
    
    return result;
}

void TextureTiler::gatherTiles()
{
    for (Triangle *tri = mesh_->meshTris_;
         (tri - mesh_->meshTris_) < mesh_->numMeshTris_;
         tri++) {
        
        findOrCreateTile( tri );
    }
    
    printf("Gather Tiles: %zu unique tiles\n", numTiles_ );
    for (size_t tileNdx = 0; tileNdx < numTiles_; tileNdx++) {
        Tile *tile = tiles_[tileNdx];
        
        char buff[100];
        sprintf( buff, "dbgtiles/tile_%c%c%c.png",
                tile->edge_[0]->edgeCode_ + 'A' - 1,
                tile->edge_[1]->edgeCode_ + 'A' - 1,
                tile->edge_[2]->edgeCode_ + 'A' - 1 );
        
        tile->img_->filename_ = strdup(buff);
    }
    
    // FIXME: delete images
}

void TextureTiler::assembleTiles()
{
    uint32_t outSz = 128;
    
    // First we need to figure out how big the output image should be
    bool tilesFit = false;

    // Try to pack them into a square image this size
    while (!tilesFit)
    {
        tilesFit = true;
        
        printf("Packing size %dx%d\n", outSz, outSz );
        uint32_t packX=0, packY=0;
        uint32_t rowY = 0;
        for (size_t tileNdx = 0; tileNdx < numTiles_; tileNdx++) {
            Tile *tile = tiles_[tileNdx];
            if ( (packX + tile->img_->width_) >= outSz ) {
                // advance to next row
                packX = 0;
                packY += rowY;
                rowY = 0;
            }
            
            tile->packX_ = packX;
            tile->packY_ = packY;
            
            packX += tile->img_->width_;
            if (tile->img_->height_ > rowY) {
                rowY = tile->img_->height_;
            }
            
            if (packY + rowY > outSz) {
                // doesn't fit, start over
                tilesFit = false;
                outSz *= 2;
                break;
            }
        }
    }
    
    // now make the output image and do the pack
    // Make an output texture
    printf("Output size: %d %d\n", outSz, outSz );
    outTexture_ = new tapnik::Image( outSz, outSz );
    outTexture_->clear( 0xff7f7f7f );
    outTexture_->filename_ = strdup(outTexFilename_);

    for (size_t tileNdx = 0; tileNdx < numTiles_; tileNdx++) {
        Tile *tile = tiles_[tileNdx];
        
        outTexture_->drawImage( tile->packX_, tile->packY_, tile->img_ );
    }
}

void TextureTiler::debugDumpTiles()
{
    for (size_t tileNdx = 0; tileNdx < numTiles_; tileNdx++) {
        Tile *tile = tiles_[tileNdx];
        
        tile->paintFromSource( sourceImage_ );
        tile->debugDrawAnnotations();
        
        tile->img_->save();
    }
}

void TextureTiler::finish()
{
    assert(outTexture_);
    outTexture_->save();
}


