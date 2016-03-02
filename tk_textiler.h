//
//  tk_textile.hpp
//  tk_textile
//
//  Created by Joel Davis on 2/27/16.
//  Copyright © 2016 Joel Davis. All rights reserved.
//

#pragma once

// FIXME: replace this with a cross-platform math library
#include <GLKit/GLKMath.h>

#include <stdint.h>

#include "tk_objfile.h"

namespace tapnik
{

// RGBA image
struct Image
{
    char *filename_ = nullptr;
    uint32_t width_ = 0;
    uint32_t height_ = 0;
    uint32_t *imgdata_ = nullptr;

    Image() {};
    Image( uint32_t width, uint32_t height );
    ~Image();

    inline void drawPixel( int32_t x, int32_t y, uint32_t color=0xffffffff );
    inline uint32_t getPixel( int32_t x, int32_t y );
    void drawLine( int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color=0xffffffff );
    void drawFatLine( int32_t x1, int32_t y1, int32_t x2, int32_t y2, uint32_t color=0xffffffff );
    void drawImage( int32_t x, int32_t y, tapnik::Image *img );
    void clear( uint32_t color=0xff000000);
    
    static Image *load( const char *filename );

    void save(); // save to filename_
};

struct EdgeInfo
{
    EdgeInfo( uint32_t code, uint32_t debugColor ) :
        edgeCode_(code), debugColor_(debugColor) {}
    
    uint32_t edgeCode_;
    uint32_t debugColor_;
    
    // Where this edge lands on src image    
    GLKVector3 srcPointA_, srcPointB_;
};
    
    
// hmm, maybe tile should be 2 triangles that share an edge?
struct Tile
{
    Tile( int32_t edgeSize, int32_t margin );
    
    // Tile triangles in pixel coords
    int32_t tileA_[2], tileB_[2], tileC_[2];
    
    EdgeInfo *edge_[3];
    tapnik::Image *img_;
    
    bool flipped_[3]; // AB, BC, CA
    
    // Xform into source image
    GLKMatrix4 xform_;
    
    // packing pos in output image
    uint32_t packX_ = 0;
    uint32_t packY_ = 0;
    
    void paintFromSource( tapnik::Image *srcImage );
    
    void debugDrawAnnotations();
};
    

struct Triangle
{
    // just reuse TK_TriangleVert since we need the same info
    TK_TriangleVert A_,B_,C_;
    
    EdgeInfo *ab_=nullptr;
    EdgeInfo *bc_=nullptr;
    EdgeInfo *ca_=nullptr;
    
    // neighbor across edge
    Triangle *nbAB_=nullptr;
    Triangle *nbBC_=nullptr;
    Triangle *nbCA_=nullptr;
    
    bool flipped_[3]; // AB, BC, CA
    
    Tile *tile_=nullptr;
    
    bool visited_;
};

#define TK_NUM_EDGE_COLORS (5)
    
struct Mesh
{
    char *filename_=nullptr;
    
    size_t numMeshTris_;
    Triangle *meshTris_;

    void buildAdjacency();
    void assignEdges( EdgeInfo *edges[TK_NUM_EDGE_COLORS] );

    size_t solvedTris_;
    void doAssign( Triangle *tri, EdgeInfo *edges[TK_NUM_EDGE_COLORS] );
    
    static Mesh *load( const char *filename );
    void save( const char *filename, int32_t outSz );
    
    ~Mesh();
};
  
    
struct TextureTiler
{
    ~TextureTiler();
    
    void doStuff( const char *outTexFilename );
    
    void assignEdges();
    void gatherTiles();
    void assembleTiles(); 
    void finish();
    
    void debugDumpTiles();
    
    void placeEdge( EdgeInfo *edge );
    
    tapnik::Tile *findOrCreateTile( tapnik::Triangle *tri );
    
    tapnik::Image *sourceImage_ = nullptr;

    char *outTexFilename_ = nullptr;
    tapnik::Image *outTexture_ = nullptr;
    
    uint32_t edgeSize_ = 160; // size of output triangle edges
    
    tapnik::Mesh *mesh_;

    size_t numTiles_=0;
    size_t tilesCapacity_=0;
    tapnik::Tile **tiles_=nullptr;
};
    
} // namespace tapnik
