/*
	Simple TGA
	Limited edition(BGRA32bit only)
 */

#ifndef SIMPLE_TGA_HEADER
#define SIMPLE_TGA_HEADER

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __GNUC__
#pragma pack(push, 1)
#define ALIGNMENT __attribute__((packed))
#else
#pragma pack(1)
#define ALIGNMENT
#endif // __GNUC__
typedef struct
{
	unsigned char identsize;
	unsigned char colourmaptype;
	unsigned char imagetype;
	short colourmapstart;
	short colourmaplength;
	unsigned char colourmapbits;
	short xstart;
	short ystart;
	short width;
	short height;
	unsigned char bits;
	unsigned char descriptor;
} TGA_HEADER;
#ifdef __GNUC__
#pragma pack(pop)
#else
#pragma pack()
#endif // __GNUC__

inline bool isLittleEndian()
{
    const int x = 1;
    return (*(char*)&x) > 0;
}

inline short swapshortEndian(short s)
{
    return ((s & 0xFF) << 8) | (s>>8);
}
void swapTGAHeaderEndian(TGA_HEADER* header)
{
    header->colourmapstart  = swapshortEndian(header->colourmapstart);
    header->colourmaplength = swapshortEndian(header->colourmaplength);
    header->xstart = swapshortEndian(header->xstart);
    header->ystart = swapshortEndian(header->ystart);
    header->width  = swapshortEndian(header->width);
    header->height = swapshortEndian(header->height);
}

bool SimpleTGALoader(const void* tgabuffer, int& w, int& h, const unsigned char** bgra)
{
	const TGA_HEADER* header = reinterpret_cast<const TGA_HEADER*>(tgabuffer);
	if (isLittleEndian()) {
	    w = header->width;
	    h = header->height;
	} else {
	    w = swapshortEndian(header->width);
	    h = swapshortEndian(header->height);
	}

	if (header->bits != 32 || header->colourmaptype != 0 || header->imagetype != 2)
		return false;
	*bgra = static_cast<const unsigned char*>(tgabuffer);
	*bgra += sizeof(TGA_HEADER);
	return true;
}

bool SimpleTGALoaderRGBA(const void* tgabuffer, int& w, int& h, unsigned char** rgba)
{
	const unsigned char* bgra;
	const bool r = SimpleTGALoader(tgabuffer, w, h, &bgra);
	if (!r)
		return false;
	*rgba = new unsigned char[4 * w * h];
	const int wxh = w * h;
	for (int i = 0; i < wxh; ++i) { // bgra -> rgba
		(*rgba)[4 * i + 0] = bgra[4 * i + 2];
		(*rgba)[4 * i + 1] = bgra[4 * i + 1];
		(*rgba)[4 * i + 2] = bgra[4 * i + 0];
		(*rgba)[4 * i + 3] = bgra[4 * i + 3];
	}
	return true;
}

int SimpleTGASaver(void** tgabuffer, int w, int h, const unsigned char* bgra)
{
	if (!tgabuffer)
		return 0;
	*tgabuffer = malloc(sizeof(TGA_HEADER) + w * h * 4);
	const TGA_HEADER header = {0,0,2, 0,0,0, 0,0,static_cast<short>(w),static_cast<short>(h),32};
	memcpy(*tgabuffer, static_cast<const void*>(&header), sizeof(TGA_HEADER));
	memcpy(static_cast<char*>(*tgabuffer)+18, bgra, w * h * 4);
	if (!isLittleEndian())
	    swapTGAHeaderEndian(static_cast<TGA_HEADER*>(*tgabuffer));
	return sizeof(TGA_HEADER) + w * h * 4;
}


int SimpleTGASaverRGBA(void** tgabuffer, int w, int h, const unsigned char* rgba)
{
	unsigned char* bgra = new unsigned char[w * h * 4];
	const int wxh = w * h;
	for (int i = 0; i <wxh; ++i) { // rgba -> bgra
		bgra[4 * i + 0] = rgba[4 * i + 2];
		bgra[4 * i + 1] = rgba[4 * i + 1];
		bgra[4 * i + 2] = rgba[4 * i + 0];
		bgra[4 * i + 3] = rgba[4 * i + 3];
	}
	int r = SimpleTGASaver(tgabuffer, w, h, bgra);
	delete [] bgra;
	return r;
}
#endif
