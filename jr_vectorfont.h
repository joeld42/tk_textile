#ifndef VECTOR_FONT_H
#define VECTOR_FONT_H

#include "tk_textiler.h"

class VectorFont
{
public:
    virtual void move( uint32_t x, uint32_t y )=0;
    virtual void pen( uint32_t color )=0 ;
    
    virtual void vprintf( tapnik::Image *img, const char *fmt,...) = 0;
    virtual void vputc( tapnik::Image *img, char c) = 0;
};


//VectorFont * createVectorFont(void);
VectorFont * createVectorFont();
void         releaseVectorFont(VectorFont *vf);

#endif
