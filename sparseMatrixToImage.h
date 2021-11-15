/*
 * Copyright Andrea Di Iorio 2021
 * This file is part of sparseMatrixToImage
 * sparseMatrixToImage is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * sparseMatrixToImage is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with sparseMatrixToImage.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SPMAT_IMG
#define SPMAT_IMG


#define PPM_HEADER_MAX_LEN 100
#define PPM_CHAN_NUM       3
#define RGB_BLACK          0xFFFFFF

#define PATTERN_DENSE      "PATTERN_DENSE"
#define MM_COO             "MM_COO"
#define DFLT_OUTPATH       TMPDIR "mat.ppm"
typedef struct{
    char* data;     //width*height*3 triples of RGB pixel
    ulong width;

    ulong height;
    char header[PPM_HEADER_MAX_LEN];
}   ppmData;

#define NNZ_PIXEL_COLOR  255
#define Z_PIXEL_COLOR    0
/*
 * convert dense matrix @mat into ppm RGB triple pixels ( in @data ) 
 * with a black dot per NZ elem 
 * assign each consecutive @step^2 matrix elements to a dot in the PPM image
 * if at least 1 nz is in the square 
 * the dot will also include unifyNearW^2 image pixel to let the dot be more visible
 */
void denseMatrixToPPM(ppmData* data,
  ulong M, ulong N, double mat[][N],uint step, ushort unifyNearW);


/*
 * map each nz eleme in @sparseMat into a pixel dots square
 * similarly as before in _denseMatrixToPPM
 */
int sparseMatrixToPPM(ppmData* data,spmat* mat,uint step, ushort unifyNearW);

#endif
