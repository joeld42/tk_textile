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

#include "jr_vectorfont.h"

// TODO list --------
// - correspondence points in placeEdge
// - graphcut between tiles
// - handle open meshes
// - make all things cmd line options
// - 

#define TK_PI (3.1415926535897932384626433832795)
#define TK_DEG2RAD (TK_PI/180.0)
#define TK_RAD2DEG (180.0/TK_PI)
#define TK_SQRT3_OVER_2 (0.86602540378443864676372)
#define TK_SGN(x) ((x<0)?-1:((x>0)?1:0))
#define TK_ABS(x) (((x)<0)?-(x):(x))
#define TK_LERP(a,b,t) (((1.0-t)*(a)) + (t*(b)))

static VectorFont *g_vectorFont = NULL;

using namespace tapnik;

// ==================================================
# pragma mark - Misc
// ==================================================

static inline GLKVector3 lerpVec3( GLKVector3 a, GLKVector3 b, float t )
{
    return GLKVector3Make( TK_LERP(a.x,b.x,t),
                           TK_LERP(a.y,b.y,t),
                           TK_LERP(a.z,b.z,t) );
}

static inline GLKVector4 lerpVec4( GLKVector4 a, GLKVector4 b, float t )
{
    return GLKVector4Make( TK_LERP(a.x,b.x,t),
                          TK_LERP(a.y,b.y,t),
                          TK_LERP(a.z,b.z,t),
                          TK_LERP(a.w,b.w,t) );
}

static inline float saturate( float x )
{
    if (x < 0.0) return 0.0;
    else if (x > 1.0) return 1.0;
    else return x;
}


static inline float smoothstep(float edge0, float edge1, float x)
{
    // Scale, bias and saturate x to 0..1 range
    x = saturate((x - edge0)/(edge1 - edge0));
    // Evaluate polynomial
    return x*x*(3 - 2*x);
}

static inline uint32_t min_u32( uint32_t a, uint32_t b ) {
    if (a < b) return a;
    else return b;
}

float randUniform()
{
    return (float)rand() / (float)RAND_MAX;
}

float randUniform( float minVal, float maxVal )
{
    return minVal + (randUniform() * (maxVal-minVal));
}

//http://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
GLKVector3 barycentric( GLKVector3 a, GLKVector3 b, GLKVector3 c,
                       GLKVector3 p )
{
    GLKVector3 result;
    
    GLKVector3 v0 = GLKVector3Subtract(b, a);
    GLKVector3 v1 = GLKVector3Subtract(c, a);
    GLKVector3 v2 = GLKVector3Subtract(p, a);
    
    float d00 = GLKVector3DotProduct( v0, v0 );
    float d01 = GLKVector3DotProduct( v0, v1 );
    float d11 = GLKVector3DotProduct( v1, v1 );
    float d20 = GLKVector3DotProduct( v2, v0 );
    float d21 = GLKVector3DotProduct( v2, v1 );
    
    float denom = d00*d11 - d01*d01;
    result.y = (d11*d20 - d01*d21) / denom;
    result.z = (d00*d21 - d01*d20) / denom;
    result.x = 1.0f - result.y - result.z;
    
    return result;
}

static inline GLKVector4 floatColor( uint32_t pixelColor )
{
    // abgr(uint32) -> rgba (float)
    return GLKVector4Make( (float)(pixelColor  & 0xff) / 255.0,
                           (float)((pixelColor >> 8) & 0xff) / 255.0,
                           (float)((pixelColor >> 16) & 0xff) / 255.0,
                           (float)((pixelColor >> 24) & 0xff) / 255.0 );
}

static inline uint32_t pixelColor( GLKVector4 floatColor )
{
    // rgba(float -> abgr(uint32)
    return ((int)(floatColor.a*0xff) << 24) |
           ((int)(floatColor.b*0xff) << 16) |
           ((int)(floatColor.g*0xff) << 8) |
            (int)(floatColor.r*0xff);
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

// for debugging use only
inline void Image::drawPixelTinted( int32_t x, int32_t y, uint32_t origColor, uint32_t tintColor, float tintAmt  )
{
    GLKVector4 origColorf = floatColor( origColor );
    GLKVector4 tintColorf = floatColor(tintColor );

    GLKVector4 colorf = lerpVec4( origColorf, tintColorf, tintAmt );
    uint32_t color = pixelColor( colorf );
    
    if ((x>=0) && (y>=0) && (x < width_) && (y < height_)) {
        imgdata_[ (y * width_) + x ] = color;
    }
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
//    return;
//    printf("drawline %d %d -> %d %d\n", x1, y1, x2, y2 );
    
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
    for (int i=-2; i <= 2; i++) {
        for( int j=-2; j <= 2; j++) {
            drawLine( x1+i, y1+j, x2+i, y2+j, color );
        }
    }
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
//        printf("Wrote output image: %s\n", filename_ );
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
    t->dbgIndex_ = (int)(mesh->numMeshTris_ - 1);
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

void dbgPrintTri( Triangle *tri )
{
    printf("Index: %d packFlip %d packRot %d\n", tri->dbgIndex_, tri->packFlip_?1:0, tri->packRot_ );
    
}

void Mesh::save( const char *filename, int32_t outMapSize )
{
    printf("Writing output mesh '%s'\n", filename );
    
    FILE *fp = fopen(filename, "wt");
    
    fprintf( fp, "# Saved from tk_textile\n");
    fprintf( fp, "usemtl testoutmap\n");
    // verts
    for (Triangle *tri = meshTris_; (tri - meshTris_) < numMeshTris_; tri++) {
        fprintf( fp, "v %f %f %f\n", tri->A_.pos[0], tri->A_.pos[1],tri->A_.pos[2]);
        fprintf( fp, "v %f %f %f\n", tri->B_.pos[0], tri->B_.pos[1],tri->B_.pos[2]);
        fprintf( fp, "v %f %f %f\n\n", tri->C_.pos[0], tri->C_.pos[1],tri->C_.pos[2]);
    }

    // normals
    for (Triangle *tri = meshTris_; (tri - meshTris_) < numMeshTris_; tri++) {
        fprintf( fp, "vn %f %f %f\n", tri->A_.nrm[0], tri->A_.nrm[1],tri->A_.nrm[2]);
        fprintf( fp, "vn %f %f %f\n", tri->B_.nrm[0], tri->B_.nrm[1],tri->B_.nrm[2]);
        fprintf( fp, "vn %f %f %f\n\n", tri->C_.nrm[0], tri->C_.nrm[1],tri->C_.nrm[2]);
    }

    // tex coords from tiles
    for (Triangle *tri = meshTris_; (tri - meshTris_) < numMeshTris_; tri++) {
        Tile *tile = tri->tile_;
        
        int pr = tri->packFlip_?tri->packRot_:(3-tri->packRot_);
        GLKVector3 tileST[3];
        if (!tri->packFlip_) {
            tileST[(pr+0)%3] = GLKVector3Make( tile->tileA_[0], tile->tileA_[1], 0.0);
            tileST[(pr+1)%3] = GLKVector3Make( tile->tileB_[0], tile->tileB_[1], 0.0);
            tileST[(pr+2)%3] = GLKVector3Make( tile->tileC_[0], tile->tileC_[1], 0.0);
//            printf ("packRot %d pr %d triA %d triB %d triC %d\n",
//                    tri->packRot_, pr,
//                    (pr+0)%3, (pr+1)%3, (pr+2)%3 );
        } else {
            tileST[(pr+0)%3] = GLKVector3Make( tile->tileA_[0], tile->tileA_[1], 0.0);
            tileST[(pr+1)%3] = GLKVector3Make( tile->tileC_[0], tile->tileC_[1], 0.0);
            tileST[(pr+2)%3] = GLKVector3Make( tile->tileB_[0], tile->tileB_[1], 0.0);
        }
        
        float stA[2], stB[2], stC[2];
        
        //        if ( tile->dbgIndex_ == 33 ) {
        if (tri->dbgIndex_ == 95) {
            //            printf("Suspicious tile on triangle %d\n", tri->dbgIndex_);
            //            tileST[1] = tileST[0];
            //            tileST[2] = tileST[0];
            
            // print out the bad tile
            printf("-------\n");
            printf("BAD STs: A %3.2f %3.2f B %3.2f %3.2f C %3.2f %3.2f\n",
                   tileST[0].x, tileST[0].y,
                   tileST[1].x, tileST[1].y,
                   tileST[2].x, tileST[2].y );
            dbgPrintTri( tri );
        }


        stA[0] = (float)(tile->packX_ + tileST[0].x ) / (float)outMapSize;
        stA[1] = 1.0 - (float)(tile->packY_ + tileST[0].y ) / (float)outMapSize;
        fprintf( fp, "vt %f %f\n", stA[0], stA[1] );
        
        stB[0] = (float)(tile->packX_ + tileST[1].x ) / (float)outMapSize;
        stB[1] = 1.0 - (float)(tile->packY_ + tileST[1].y ) / (float)outMapSize;
        fprintf( fp, "vt %f %f\n", stB[0], stB[1] );

        stC[0] = (float)(tile->packX_ + + tileST[2].x ) / (float)outMapSize;
        stC[1] = 1.0 - (float)(tile->packY_ + tileST[2].y ) / (float)outMapSize;
        fprintf( fp, "vt %f %f\n\n", stC[0], stC[1] );
        
        
    }

    // And finally triangles
    size_t ndx = 1;
    for (Triangle *tri = meshTris_; (tri - meshTris_) < numMeshTris_; tri++) {
        fprintf( fp, "f %zu/%zu/%zu %zu/%zu/%zu %zu/%zu/%zu\n",
                ndx, ndx, ndx,
                ndx+1, ndx+1, ndx+1,
                ndx+2, ndx+2, ndx+2 );
        ndx += 3;
    }
    fclose( fp );
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

// should be factorial of TK_MAX_EDGE_COLORS + 1 (the extra is just padding for a copy)
#define TK_NUM_EDGE_PERMS (121)
static int edgePerms[TK_NUM_EDGE_PERMS][TK_MAX_EDGE_COLORS];
static int nEdgePerms = 0;

static void assignBackEdge( Triangle *tri, Triangle *nbr, EdgeInfo *edge, bool flipped )
{
    if (nbr->nbAB_ == tri) {
        nbr->ab_ = edge;
        nbr->flipped_[0] = flipped;
    } else if (nbr->nbBC_ == tri) {
        nbr->bc_ = edge;
        nbr->flipped_[1] = flipped;
    } else if (nbr->nbCA_ == tri) {
        nbr->ca_ = edge;
        nbr->flipped_[2] = flipped;
    } else {
        printf("WARN: couldn't find back edge for nbr.\n");
    }
}

static inline bool checkNeighbor( Triangle *tri, int numEdgeColors ) {
    int count[TK_MAX_EDGE_COLORS] = {0};
    
    if (tri->ab_) count[tri->ab_->edgeCode_]++;
    if (tri->bc_) count[tri->bc_->edgeCode_]++;
    if (tri->ca_) count[tri->ca_->edgeCode_]++;
    
    for (int i=0; i < numEdgeColors; i++) {
        if (count[i] > 1) {
//            printf("Edge %d used %d times...\n", i, count[i]);
            return false;
        }
    }
    return true;
}

void Mesh::doAssign( Triangle *tri, EdgeInfo *edges[TK_MAX_EDGE_COLORS], int numEdgeColors )
{
    // Pick a random order to assign in
    int p = (int)randUniform(0, TK_NUM_EDGE_PERMS-1);
    
    for (int xx = 0; xx < numEdgeColors; xx++) {
        for (int yy = 0; yy < numEdgeColors; yy++) {
            if (yy==xx) continue;
            for (int zz = 0; zz < numEdgeColors; zz++) {
                if ((zz==xx) || (zz==yy)) continue;
                
                // get permuted edge colors
                int x = edgePerms[p][xx];
                int y = edgePerms[p][yy];
                int z = edgePerms[p][zz];
                
                // now we are trying to assign edges X, Y, Z to tri.
                // first, make sure any pre-assigned edges don't conflict.
                if ((tri->ab_) && (tri->ab_->edgeCode_ != x)) continue;
                if ((tri->bc_) && (tri->bc_->edgeCode_ != y)) continue;
                if ((tri->ca_) && (tri->ca_->edgeCode_ != z)) continue;
                
                bool assignedAB = false;
                bool assignedBC = false;
                bool assignedCA = false;
                
                // FIXME: handle triangles with open edges (null neighbors)
                
                // No previous assignment conflicts, go ahead and assign our edges
                bool badAssign = false;
                if (!tri->ab_) {
                    assignedAB = true;
                    tri->ab_ = edges[x];
                    assignBackEdge( tri, tri->nbAB_, tri->ab_, !tri->flipped_[0] );
                    
                    if (!checkNeighbor(tri->nbAB_, numEdgeColors)) badAssign = true;
                    
                }
                if ((!tri->bc_) && (!badAssign)) {
                    assignedBC = true;
                    tri->bc_ = edges[y];
                    assignBackEdge( tri, tri->nbBC_, tri->bc_, !tri->flipped_[1] );
                    
                    if (!checkNeighbor(tri->nbBC_, numEdgeColors)) badAssign = true;
                }
                if ((!tri->ca_) && (!badAssign)) {
                    assignedCA = true;
                    tri->ca_ = edges[z];
                    assignBackEdge( tri, tri->nbCA_, tri->ca_, !tri->flipped_[2] );
                    
                    if (!checkNeighbor(tri->nbCA_, numEdgeColors)) badAssign = true;
                }
                
                if (!badAssign)
                {
                    // ok, we've assigned all our edges, call our neighbors to assign them
                    tri->visited_ = true;
                    solvedTris_++;
                    
                    static int logcount=0;
                    if (logcount++>=1000) {
                        printf("solved: %zu/%zu\n", solvedTris_, numMeshTris_);
                        logcount = 0;
                    }
                    
                    // If we've solved all but the last tri, we've really solved the whole thing because we
                    // know it's a closed mesh and that last triangle will be colored by it's neighbors
                    // FIXME: generalize this by checking solved by counting unique edges not triangles.
                    if (solvedTris_==numMeshTris_-1) {
                        solvedTris_ = numMeshTris_;
                    }
                    
                    // Did we find a solution?
                    if (solvedTris_ == numMeshTris_) return;
                    
                    if (!tri->nbAB_->visited_) doAssign( tri->nbAB_, edges, numEdgeColors );
                    if (solvedTris_ == numMeshTris_) return;
                    
                    if (!tri->nbBC_->visited_) doAssign( tri->nbBC_, edges, numEdgeColors );
                    if (solvedTris_ == numMeshTris_) return;
                    
                    if (!tri->nbCA_->visited_) doAssign( tri->nbCA_, edges, numEdgeColors );
                    if (solvedTris_ == numMeshTris_) return;
                    
                    // Couldn't find a solution, back out and try some more
                    solvedTris_--;
                    tri->visited_ = false;
                }
                
                // back out any assigments that we did
                if (assignedAB) {
                    tri->ab_ = nullptr;
                    assignBackEdge( tri, tri->nbAB_, nullptr, false );
                }

                if (assignedBC) {
                    tri->bc_ = nullptr;
                    assignBackEdge( tri, tri->nbBC_, nullptr, false );
                }

                if (assignedCA) {
                    tri->ca_ = nullptr;
                    assignBackEdge( tri, tri->nbCA_, nullptr, false );
                }
            }
        }
    }
}


void makePermTable( int ndx, bool *used, int numEdgeColors )
{
    if (ndx==numEdgeColors) {
        // found a permutation, yay
//        printf("PERM %d: ", nEdgePerms );
        for (int i = 0; i < ndx; i++) {
//            printf("%d ", edgePerms[nEdgePerms][i] );
            edgePerms[nEdgePerms+1][i] = edgePerms[nEdgePerms][i];
        }
//        printf("\n");
        nEdgePerms++;
    } else {
        for (int i = 0; i < numEdgeColors; i++) {
            if (!used[i]) {
                used[i] = true;
                edgePerms[nEdgePerms][ndx] = i;
                makePermTable( ndx+1, used, numEdgeColors );
                used[i] = false;
            }
        }
    }
}

void Mesh::assignEdges( EdgeInfo *edges[TK_MAX_EDGE_COLORS], int numEdgeColors )
{
    printf("Assign edges\n");

    // generate permutation table
    nEdgePerms = 0;
    bool used[TK_MAX_EDGE_COLORS] = {0};
    makePermTable( 0, used, numEdgeColors );

    for (Triangle *tri = meshTris_; (tri - meshTris_) < numMeshTris_; tri++) {
        tri->visited_ = false;
    }
    
    solvedTris_ = 0;
    doAssign( meshTris_, edges, numEdgeColors );
    
    if (solvedTris_ == numMeshTris_) {
        printf("Yay found a solution...\n");
        
    } else {
        printf("No edge coloring found. :(\n");
    }
    
    for (Triangle *tri = meshTris_; (tri - meshTris_) < numMeshTris_; tri++) {
        // for now, just pick at random
        if (!tri->ab_) {
            tri->ab_ = edges[ rand() % 3 ];
            assignBackEdge( tri, tri->nbAB_, tri->ab_, !tri->flipped_[0] );
        }
        
        if (!tri->bc_) {
            tri->bc_ = edges[ rand() % 3 ];
            assignBackEdge( tri, tri->nbBC_, tri->bc_,  !tri->flipped_[1] );
        }
        
        if (!tri->ca_) {
            tri->ca_ = edges[ rand() % 3 ];
            assignBackEdge( tri, tri->nbCA_, tri->ca_,  !tri->flipped_[2] );
        }
    }
    
    printf("Assign edges done...\n");
}

// ==================================================
# pragma mark - Tile
// ==================================================
Tile::Tile()
{
    img_ = nullptr;
}

void Tile::makeImage( int32_t edgeSize, int32_t margin ) {
    
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
    
    g_vectorFont->pen( 0xffffffff );
    g_vectorFont->move( tileA_[0]-4, tileA_[1]+8 );
    g_vectorFont->vprintf( img_, "A" );
    
    g_vectorFont->move( tileB_[0]-21, tileB_[1] - 20 );
    g_vectorFont->vprintf( img_, "B");
    
    g_vectorFont->move( tileC_[0]+12, tileC_[1] - 20 );
    g_vectorFont->vprintf( img_, "C %s", dbgIndexStr );

//    g_vectorFont->move( (img_->width_/2) - 20, img_->height_/2 );
//    g_vectorFont->vprintf( img_, "%s%s%s", flipped_[0]?"Y":"N", flipped_[1]?"Y":"N", flipped_[2]?"Y":"N" );
//
//    g_vectorFont->move( (img_->width_/2) - 20, img_->height_/2 + 20);
//    g_vectorFont->vprintf( img_, "FARE" );

}

void combineBlend( GLKVector3 ta, GLKVector3 tb, GLKVector3 tc,
                  Image *img, Image *overImg, bool leftCorner )
{
    float blendSharpness = 0.9; // higher value = more sharper edge
    // blend
    for (int j=0; j < img->height_; j++) {
        for (int i=0; i < img->width_; i++) {
            
            GLKVector3 p = GLKVector3Make( (float)i, (float)j, 0.0 );
            GLKVector3 b = barycentric(ta, tb, tc, p );
            
            float ov = 0.1;
            if ((b.x >=-ov) && (b.x <=1.0+ov) &&
                (b.y >=-ov) && (b.y <=1.0+ov) &&
                (b.z >=-ov) && (b.z <=1.0+ov) )
            {
                b.x = saturate(b.x);
                b.y = saturate(b.y);
                b.z = saturate(b.z);
                float aa = smoothstep( 0.0, leftCorner?b.x:b.y, b.z );
                float bb = smoothstep( 0.0, leftCorner?b.x:b.y, leftCorner?b.y:b.x );
                float blendVal = pow(aa*bb, blendSharpness);
                
                GLKVector4 baseColorf = floatColor( img->getPixel(i,j) );
                GLKVector4 tileColorf = floatColor( overImg->getPixel(i, j));
                GLKVector4 blendColorf = lerpVec4( baseColorf, tileColorf, blendVal );
                blendColorf.a = 1.0; // fixme
                uint32_t blendColor = pixelColor(blendColorf );
                
                img->drawPixel(i, j, blendColor );
            }
            
        }
    }
    
}

void Tile::paintFromSource(Image *srcImage)
{
    // Paint first edge onto the tile
//    int ndx = 0; // dbg crap
//    int targ = 3;
//    if (edge_[1]->edgeCode_ == targ) ndx = 1;
//    else if (edge_[2]->edgeCode_ == targ) ndx = 2;
    paintFromSourceEdge( img_, srcImage, 0 );
    
    // Paint the next edge onto a temp image
    Image *tmpImg = new Image( img_->width_, img_->height_ );
    paintFromSourceEdge( tmpImg, srcImage, 1 );
    
    GLKVector3 ta = GLKVector3Make( tileA_[0], tileA_[1], 0.0 );
    GLKVector3 tb = GLKVector3Make( tileB_[0], tileB_[1], 0.0 );
    GLKVector3 tc = GLKVector3Make( tileC_[0], tileC_[1], 0.0 );

    combineBlend(ta, tb, tc, img_, tmpImg, true );
    

    // Paint the last edge onto a temp image
    paintFromSourceEdge( tmpImg, srcImage, 2 );

    combineBlend(ta, tb, tc, img_, tmpImg, false );
    
    delete tmpImg;
}

void Tile::paintFromSourceEdge(Image *destImage, Image *srcImage, int edgeIndex )
{
    // Get requested edge
    EdgeInfo *edge;
    GLKVector3 a, b;
    if (edgeIndex==0) {
        //edge 0 is AB
        edge = edge_[0];
        if (!flipped_[0]) {
            a = GLKVector3Make( tileA_[0], tileA_[1], 0.0 );
            b = GLKVector3Make( tileB_[0], tileB_[1], 0.0 );
        } else {
            b = GLKVector3Make( tileA_[0], tileA_[1], 0.0 );
            a = GLKVector3Make( tileB_[0], tileB_[1], 0.0 );
        }
    } else if (edgeIndex==1) {
        //edge 1 is BC
        edge = edge_[1];
        if (!flipped_[1]) {
            a = GLKVector3Make( tileB_[0], tileB_[1], 0.0 );
            b = GLKVector3Make( tileC_[0], tileC_[1], 0.0 );
        } else {
            b = GLKVector3Make( tileA_[0], tileA_[1], 0.0 );
            a = GLKVector3Make( tileB_[0], tileB_[1], 0.0 );
        }
    } else { // edgeIndex==2
        //edge 2 is CA
        edge = edge_[2];
        if (!flipped_[2]) {
            a = GLKVector3Make( tileC_[0], tileC_[1], 0.0 );
            b = GLKVector3Make( tileA_[0], tileA_[1], 0.0 );
        } else {
            b = GLKVector3Make( tileC_[0], tileC_[1], 0.0 );
            a = GLKVector3Make( tileA_[0], tileA_[1], 0.0 );
        }
    }

    
    GLKVector3 destA = edge->srcPointB_;
    GLKVector3 destB = edge->srcPointA_;
    GLKVector3 srcEdgeDir = GLKVector3Subtract( b, a );
    GLKVector3 destEdgeDir = GLKVector3Subtract( destB, destA );
    
    float angSrc =atan2f( srcEdgeDir.y, srcEdgeDir.x );
    float angDest = atan2f( destEdgeDir.y, destEdgeDir.x );
    
    
    GLKVector3 translate = GLKVector3Subtract( destA, a );
    float lengthAB = GLKVector3Length( GLKVector3Subtract( a, b ) );
    float lengthDestAB = GLKVector3Length( GLKVector3Subtract( destA, destB ) );
    float scale = lengthAB / lengthDestAB;
    float angle = angDest - angSrc;
    
    // scale should be ~= 1.0 since we generated the target line the same size
//    printf("lengthA %f lengthDest %f Scale is %f\n", lengthAB, lengthDestAB, scale );
    
    GLKMatrix4 m3 = GLKMatrix4MakeTranslation( a.x, a.y, a.z );
    GLKMatrix4 m2 = GLKMatrix4MakeRotation( angle, 0.0, 0.0, 1.0 );
    GLKMatrix4 m1 = GLKMatrix4MakeTranslation( -a.x, -a.y, -a.z );
    
    GLKMatrix4 xform1 = GLKMatrix4Multiply(m3, GLKMatrix4Multiply(m2, m1));
    GLKMatrix4 xform2 = GLKMatrix4MakeTranslation( translate.x, translate.y, translate.z );
    
//    xform_ = xform2;
    xform_ = GLKMatrix4Multiply( xform2, xform1 );
    
    // DBG: draw transformed triangle into source image
    GLKVector3 ta = GLKVector3Make( tileA_[0], tileA_[1], 0.0 );
    GLKVector3 tb = GLKVector3Make( tileB_[0], tileB_[1], 0.0 );
    GLKVector3 tc = GLKVector3Make( tileC_[0], tileC_[1], 0.0 );
//    GLKVector3 aa = GLKMatrix4MultiplyVector3WithTranslation( xform_, GLKVector3Make( tileA_[0], tileA_[1], 0.0 ));
//    GLKVector3 bb = GLKMatrix4MultiplyVector3WithTranslation( xform_, GLKVector3Make( tileB_[0], tileB_[1], 0.0 ));
//    GLKVector3 cc = GLKMatrix4MultiplyVector3WithTranslation( xform_, GLKVector3Make( tileC_[0], tileC_[1], 0.0 ));
    
//    srcImage->drawLine( aa.x, aa.y, bb.x, bb.y );
//    srcImage->drawLine( bb.x, bb.y, cc.x, cc.y );
//    srcImage->drawLine( cc.x, cc.y, aa.x, aa.y );
    
    // Copy the pixels from source into tile
    for (int j=0; j < img_->height_; j++) {
        for (int i=0; i < img_->width_; i++) {
            
            GLKVector3 p = GLKVector3Make( (float)i, (float)j, 0.0 );
            GLKVector3 b = barycentric(ta, tb, tc, p );
            
            float ov = 0.1;
            if ((b.x >=-ov) && (b.x <=1.0+ov) &&
                (b.y >=-ov) && (b.y <=1.0+ov) &&
                (b.z >=-ov) && (b.z <=1.0+ov) )
            {
                GLKVector4 samplePos = GLKVector4Make( (float)i, (float)j, 0.0, 1.0 );
                samplePos = GLKMatrix4MultiplyVector4( xform_, samplePos);
                
                // TODO: (maybe) fractional lookup and interpolate
                uint32_t sampleVal = srcImage->getPixel( (int32_t)samplePos.x, (int32_t)samplePos.y );
#if 0
                destImage->drawPixel(i, j, sampleVal );
#else
                destImage->drawPixelTinted(i, j, sampleVal, edge->debugColor_, 0.5 );
#endif
            }
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

void TextureTiler::placeEdge( EdgeInfo *edge )
{
    float edgeSz = (float)edgeSize_;
    float halfEdgeSz = edgeSz / 2.0;
    float randAngle = randUniform( 0.0, 360.0 ) * TK_DEG2RAD;
    
    GLKVector3 center = GLKVector3Make( randUniform( edgeSz, sourceImage_->width_ - edgeSz) ,
                                        randUniform( edgeSz, sourceImage_->height_ - edgeSz)  , 0.0 );
    
    GLKVector3 v = GLKVector3MultiplyScalar( GLKVector3Make( cosf(randAngle), sinf(randAngle), 0.0 ), halfEdgeSz );
    
    edge->srcPointA_ = GLKVector3Add( center, v );
    edge->srcPointB_ = GLKVector3Subtract( center, v );
    
    printf("Edge %d -- A %3.2f %3.2f %3.2f B %3.2f %3.2f %3.2f\n",
           edge->edgeCode_,
           edge->srcPointA_.x, edge->srcPointA_.y, edge->srcPointA_.z,
           edge->srcPointB_.x, edge->srcPointB_.y, edge->srcPointB_.z );
}

// Does stuff.
void TextureTiler::doStuff( const char *outTexFilename )
{
    g_vectorFont = createVectorFont();

//    g_vectorFont->pen( 0xffff00ff );
//    g_vectorFont->move( 100, 100 );
//    g_vectorFont->vprintf( sourceImage_, "Hello world" );
    
    // Make edge colors
    EdgeInfo *edges[TK_MAX_EDGE_COLORS];
    edges[0] = new EdgeInfo( 0, 0xffff0000 );
    edges[1] = new EdgeInfo( 1, 0xff00ff00 );
    edges[2] = new EdgeInfo( 2, 0xff0000ff );
    edges[3] = new EdgeInfo( 3, 0xff008efc );
    edges[4] = new EdgeInfo( 4, 0xff7f007f );

    // build mesh
    mesh_->buildAdjacency();
    mesh_->assignEdges( edges, numEdgeColors_ );
    
//    // Save result
//    FILE *fp = fopen("assign.txt", "wt");
//    for (int i=0; i < numTiles_; i++) {
//        
//        fprintf(fp, "%d: %d %d %d %s%s%s\n",
//                i,
//    }
    
    // make tiles
    gatherTiles();

    // calc edgeSize
    int bestSize = -1;
    int bestCols = 0;
    float bestDiff = 0;
    for (int rowCount=1; rowCount <= numTiles_; rowCount++) {
        
        // See what the edge size would be for this rowCount
        // TODO: flip alternate triangles for tighter pack
        float szX = 0;
        float rowX = 0;
        float szY = 0;
        int r = 0;
        int numCols = 1;
        for (int i=0; i < numTiles_; i++) {
            rowX += 1;
            r++;
            if (r==rowCount) {
                r = 0;
                numCols++;
                rowX = 0;
                szY += TK_SQRT3_OVER_2;
            }
            if (rowX > szX) szX = rowX;
        }
        float diff = fabs(szX - szY);
        if ((bestSize<0)||(diff<bestDiff)) {
            printf("At rowcount %d sz is %f x %f\n", rowCount, szX, szY );
            printf("rowCount %d Diff %f bestDiff %f\n", rowCount, diff, bestDiff );
            bestDiff = diff;
            bestSize = rowCount;
            bestCols = numCols;
        }
    }
    
    int edgeSize1 = (int)((outputSize_ - 8.0) / (float)bestSize) - 4.0; // required size from width
    int edgeSize2 = (int)( floorf( (outputSize_-8.0) / (float)bestCols) * (1.0/TK_SQRT3_OVER_2)) - 8.0; // required size from height

    float res1 = edgeSize1 * bestSize;
    float res2 = (edgeSize2*TK_SQRT3_OVER_2) * bestCols;
    printf("res1: %f\n", res1);
    printf("res2: %f\n", res2);
    
    edgeSize_ = min_u32( edgeSize1, edgeSize2 );
    printf("edgeSize1 %d edgeSize2 %d\n", edgeSize1, edgeSize2 );
    printf("Best Row Size: %d bestColSize: %d, edgeSize %d\n", bestSize, bestCols, edgeSize_);
    
    // Generate images
    for (size_t tileNdx = 0; tileNdx < numTiles_; tileNdx++) {
        Tile *tile = tiles_[tileNdx];
        
        tile->makeImage( edgeSize_, 4 );
        
        char buff[100];
        sprintf( buff, "dbgtiles/tile_%c%c%c%c%c%c.png",
                tile->flipped_[0]?'n':'f', tile->edge_[0]->edgeCode_ + 'A',
                tile->flipped_[1]?'n':'f', tile->edge_[1]->edgeCode_ + 'A',
                tile->flipped_[2]?'n':'f', tile->edge_[2]->edgeCode_ + 'A' );
        
        tile->img_->filename_ = strdup(buff);
    }
    
    // DBG
//    assembleTiles(bestSize, 4 );
//    exit(1);

    
    // FIXME: delete images

    
    // initial placement of edges in source
    for (int i=0; i < numEdgeColors_; i++) {
        placeEdge( edges[i] );
        
        sourceImage_->drawFatLine( edges[i]->srcPointA_.x, edges[i]->srcPointA_.y,
                                   edges[i]->srcPointB_.x, edges[i]->srcPointB_.y,
                                  edges[i]->debugColor_ );
    }
    
    
    // paint tiles
    debugDumpTiles();
    
    assembleTiles( bestSize, 4 );

    // TODO: make filename and stuff settable
    // NOTE: this assumes outTex is square
    mesh_->save("outmesh.obj", outTexture_->width_ );
    
    finish();
}

Tile *TextureTiler::findOrCreateTile( Triangle *tri )
{
    Tile *result = nullptr;
    for (int tileRot = 0; tileRot < 3; tileRot++ )
    {
        // DBG
//        if (tileRot!=0) continue;
        
        for (size_t tileNdx = 0; tileNdx < numTiles_; tileNdx++) {
            Tile *tile = tiles_[tileNdx];

            if ( (tri->ab_ == tile->edge_[(tileRot+0)%3]) && (tri->flipped_[0] == tile->flipped_[(tileRot+0)%3]) &&
                 (tri->bc_ == tile->edge_[(tileRot+1)%3]) && (tri->flipped_[1] == tile->flipped_[(tileRot+1)%3]) &&
                 (tri->ca_ == tile->edge_[(tileRot+2)%3]) && (tri->flipped_[2] == tile->flipped_[(tileRot+2)%3])) {
                result = tile;
                tri->packFlip_ = false;
                tri->packRot_ = tileRot;
                break;
            }

            // Confusing: because packFlip reverses triangle direction, we want to match the same flip flags
            // to get the opposite
            if ( (tri->ca_ == tile->edge_[(tileRot+0)%3]) && (tri->flipped_[2] == !tile->flipped_[(tileRot+0)%3]) &&
                (tri->bc_ == tile->edge_[(tileRot+1)%3]) && (tri->flipped_[1] == !tile->flipped_[(tileRot+1)%3]) &&
                (tri->ab_ == tile->edge_[(tileRot+2)%3]) && (tri->flipped_[0] == !tile->flipped_[(tileRot+2)%3])) {
                result = tile;
                tri->packFlip_ = true;
                tri->packRot_ = tileRot;
                break;
            }

        }
    }
    
    // Didn't find one, create a new tile
    if (!result) {
        tri->packFlip_ =false;
        
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
        
        result = new Tile();
        result->edge_[0] = tri->ab_;
        result->edge_[1] = tri->bc_;
        result->edge_[2] = tri->ca_;
        
        sprintf( result->dbgIndexStr, "T%zu: %d", numTiles_, tri->dbgIndex_);
        result->dbgIndex_ = (int)numTiles_;
        
        for (int i=0; i < 3; i++) {
            result->flipped_[i] = tri->flipped_[i];
        }
        
        tiles_[numTiles_++] = result;
    }
    
    return result;
}

void TextureTiler::gatherTiles()
{
    for (Triangle *tri = mesh_->meshTris_;
         (tri - mesh_->meshTris_) < mesh_->numMeshTris_;
         tri++) {
        
        tri->tile_ = findOrCreateTile( tri );
    }
}

void TextureTiler::assembleTiles( int rowCount, int margin )
{
    int marg = 4;
    uint32_t packX=marg;
    uint32_t packY=marg;
    uint32_t rowY = 0;
    int r = 0;
    for (size_t tileNdx = 0; tileNdx < numTiles_; tileNdx++) {
        Tile *tile = tiles_[tileNdx];
        r++;
        if ( r > rowCount ) {
            // advance to next row
            packX = marg;
            packY += rowY;
            rowY = 0;
            r=1;
        }
        tile->packX_ = packX;
        tile->packY_ = packY;
        
        packX += tile->img_->width_;
        if (tile->img_->height_ > rowY) {
            rowY = tile->img_->height_;
        }
    }
    
    // now make the output image and do the pack
    // Make an output texture
    printf("Output size: %d %d\n", outputSize_, outputSize_);
    outTexture_ = new tapnik::Image( outputSize_, outputSize_ );
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
    
    // save the annotated source image
    sourceImage_->filename_ = strdup( "dbg_source.png" );
    sourceImage_->save();
}

void TextureTiler::finish()
{
    assert(outTexture_);
    outTexture_->save();
}


