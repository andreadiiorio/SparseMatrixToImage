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
 *******************************************************************************
 * dev Andrea Di Iorio
 * draw a [sparse] matrix as an image with a black square dot for each nonzero elem
 * scaling the size of the matrix, "oring" the nz elem in a square of the original matrix into 
 * the smaller destination pixel grid.
 * Highlight nz elem into square of pixel of custom size, that will incorporate near elements
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>

#include "macros.h"
#include "sparseMatrix.h"
#include "sparseMatrixToImage.h"
#include "utils.h"
#include "parser.h"

//inline export 
int BISECT_ARRAY(ulong target, ulong* arr, ulong len);
int IS_NNZ(spmat* smat,ulong i,ulong j);

///aux
double* parseDenseMatrixPatternFromFile(char* fpath,ulong* M,ulong* N){
    ERRPRINT("TODO IMPLEMENT");
    return NULL;
}
//set @i,@j in the scaled pixel grid into raw PPM -> RGB data as a black dot
static inline void _setBlackDotRGB(ulong w,ulong h, uchar rawData[][w],long i,long j){
    for ( int jj=j; jj<j+PPM_CHAN_NUM; jj++ )   rawData[i][jj] = (uchar) NNZ_PIXEL_COLOR;
    /*if  ( i<0 || i> (int)h )        return;
    for ( int jj=j; jj<j+PPM_CHAN_NUM; jj++ ){
        //skip outOfBorder enlarged dots
        if ( jj<0 || jj > (int) w)    continue; 
        //set RGB black -> no intensity in every channel
        rawData[i][jj] = NNZ_PIXEL_COLOR;
    } */
}

/*
 * set a black dot pixel in MATRIX coordinate @i,@j , that will be enlarged
 * to be more visible to a square of @unifyNearW^2 near pixels 
 * (discarding previous 0 -> white dots, 
 *  NB never setted a white dot, memset at the start)
 *  TODO FIX unifyNearW!!
 */ 
static void setBlackNZPixel(ppmData* data,long i,long j,ushort unifyNearW){
    //set the first dot
    _setBlackDotRGB(data->width*PPM_CHAN_NUM,data->height, 
      (uchar (*)[data->width*PPM_CHAN_NUM]) data->data ,i , j*PPM_CHAN_NUM);
    //set the enlarged dots
    /*for (short ww,w=0;w<unifyNearW;w++){      TODO FIX unifyNearW
        for (short zz,z=0;z<unifyNearW;z++){
            //make the highlight unify square centered in (@i,@j)
            ww = INT_DIV_CEIL(w,2); 
            zz = INT_DIV_CEIL(z,2); 
            if (!(w % 2))   ww *= -1;
            if (!(z % 2))   zz *= -1;

            _setBlackDotRGB(data->width*PPM_CHAN_NUM,data->height, 
              (uchar (*)[data->width*PPM_CHAN_NUM]) data->data, 
              i+ww,PPM_CHAN_NUM*(j+zz));
        }
    }*/
}

/////TOPLEVEL IMAGE CONVERT FUNCTIONS
void denseMatrixToPPM(ppmData* data,
  ulong M, ulong N, double mat[][N],uint step, ushort unifyNearW){
    char nz;    //flag to mark a founded nz
    for (ulong i=0;i<M;i+=step){
        for (ulong j=0;j<N;j+=step){
            //i,j point to the first: top,left element in the search square
            nz = 0;
            for (uint w=0; w<step && !nz; w++){
                for (uint z=0; z<step && !nz; z++){
                    if (mat[i+w][j+z]){
                        nz = (!0);
                        setBlackNZPixel(data,i/step,j/step,unifyNearW);
                        break;
                    }
                }
            }
        }
    }
}
int sparseMatrixToPPM(ppmData* data,spmat* mat,
    uint step, ushort unifyNearW){
    //flag matrix to set NZ element in coarsened pixel grid already colored
    char* auxPixGridFlag=calloc(data->width*data->height,sizeof(auxPixGridFlag)); 
    if (!auxPixGridFlag){
        ERRPRINT("sparseMatrixToPPM: auxPixGridFlag calloc errd\n");
        return EXIT_FAILURE;
    }
    #pragma omp parallel for schedule(static)
    for (ulong r=0;  r<mat->M; r++){
        for (ulong j=mat->IRP[r],c,r_pix,c_pix; j<mat->IRP[r+1]; j++){
            c = mat->JA[j];
            //map NZ element column into a pixel grid col
            r_pix = r/step;
            c_pix = c/step;
            //already colored this element in pixel coarsened grid
            if ( auxPixGridFlag[ IDX2D(r_pix,c_pix,data->width) ] )  continue;
            auxPixGridFlag[ IDX2D(r_pix,c_pix,data->width) ] = !0;
            setBlackNZPixel(data,r_pix,c_pix,unifyNearW);
        }
    }
    VERBOSE{
        ERRPRINT("sparseMatrixToPPM OVER");
        fflush(NULL);
    }
    _free:
    free(auxPixGridFlag);
    
    return EXIT_SUCCESS;
}

#ifdef TEST
//check if MXN dense matrix @mat has a black pixel corresponding for each nonzero element
//TODO add support for step, unifyNearW, NOW CONSIDERED AS dflt... ->1a1 mat
static int checkDenseMatrixToPPM(ulong M, ulong N,
  unsigned char rawData[][3*N],double mat[][N],ushort step, ushort unifyNearW){
    for (ulong i=0; i<M; i++){ 
        for (ulong j=0; j<N; j++){
            if (mat[i][j] && (rawData[i][3*j+0] != NNZ_PIXEL_COLOR || 
              rawData[i][3*j+1] != NNZ_PIXEL_COLOR || rawData[i][3*j+2] != NNZ_PIXEL_COLOR)){
                fprintf(stderr,"not matching NNZ at (%lu,%lu)\n",i,j);
                return EXIT_FAILURE;
            }
            else if (!mat[i][j] && (rawData[i][3*j+0] != Z_PIXEL_COLOR || 
              rawData[i][3*j+1] != Z_PIXEL_COLOR || rawData[i][3*j+2] != Z_PIXEL_COLOR)){
                fprintf(stderr,"not matching Z at (%lu,%lu)\n",i,j);
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}
//check if MXN sparse matrix @mat has a colored pixel corresponding for each nonzero element
//TODO add support for step, unifyNearW, NOW CONSIDERED AS dflt... ->1a1 mat
static int checkSparseMatrixToPPM(ulong M, ulong N,
  unsigned char rawData[][3*N],spmat* mat,ushort step, ushort unifyNearW){
    double val;
    //check NZ element in mat ARE colored in the img
    for (ulong r=0,c; r<M; r++){ 
        for (ulong j=mat->IRP[r]; j<mat->IRP[r+1]; j++){
            c   = mat->JA[j];
            val = mat->AS[j];
            if (val && (rawData[r][3*c+0] != NNZ_PIXEL_COLOR || 
              rawData[r][3*c+1] != NNZ_PIXEL_COLOR || rawData[r][3*c+2] != NNZ_PIXEL_COLOR)){
                fprintf(stderr,"not matching NNZ at (%lu,%lu)\n",r,c);
                return EXIT_FAILURE;
            }
        }
    }
    //check if other pixels are NOT colored
    for (ulong r=0; r<M; r++){ 
        for (ulong c=0,j=mat->IRP[r],cc=mat->JA[j]; c<N; c++){
            //keep checking with NZ columns to only check for Z elements
            if (c > cc && j < mat->IRP[r+1])    cc = mat->JA[++j];
            if (c == cc)    continue;   //skip NZ elements, checked above
            //ASSERT ZERO ELEMENT
            if ((rawData[r][3*c+0] != Z_PIXEL_COLOR || 
              rawData[r][3*c+1] != Z_PIXEL_COLOR || rawData[r][3*c+2] != Z_PIXEL_COLOR)){
                fprintf(stderr,"not matching Z at (%lu,%lu)\n",r,c);
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}
#endif
#ifdef MAIN_SPMAT_IMG
#define USAGE "usage: <inputMatrixFile[.gz], "PATTERN_DENSE" || "MM_COO" (matrix format) >," \
         " [MAX_OUT_SIZE,MAX_OUT_RESOLUTION_SIDE, collapseSquareW=1, outPPMFPath="DFLT_OUTPATH", unifyNearW=0]\n"
/*
 * build this file as a stand alone executable, that get the input matrix from a serialized file
 * along with the parameters from argv defining MAT2PPM_STANDALONE
 * otherwise just export the next function
 */

#include <limits.h>
//inline export 
void freeSpmatInternal(spmat*);
void freeSpmat(spmat*);
int allocSpMatrixInternal(ulong rows, ulong cols, spmat* mat);

int main(int argc,char** argv){
    void* map = NULL;
    int out=EXIT_FAILURE,outFd=0,mode=S_IRWXU;
    if (argc < 3 )    {
        ERRPRINT(USAGE);
        return EXIT_FAILURE;
    }
    char* trgtMatrix = TMP_EXTRACTED_MARTIX;
    if (extractInTmpFS(argv[1],TMP_EXTRACTED_MARTIX) < 0)   trgtMatrix = argv[1];
    double* denseMat=NULL;
    spmat*  sparseMat=NULL;
    ppmData* data=NULL;
    ulong M,N;
    if (!strncmp(argv[2],PATTERN_DENSE,strlen(PATTERN_DENSE))){
        if(!(denseMat = parseDenseMatrixPatternFromFile(argv[1],&M,&N))){
            ERRPRINT("dense matrix pattern parse failed\n");
            return EXIT_FAILURE;
        }
    } else if (!strncmp(argv[2],MM_COO,strlen(MM_COO))){
        if (!(sparseMat = MMtoCSR(trgtMatrix))){
            ERRPRINT("sparse matrix MM coordinate parsing failed\n");
            return EXIT_FAILURE;
        }
        M = sparseMat -> M;
        N = sparseMat -> N;
    }else {ERRPRINT("INVALID IN MATRIX FORMAT!\n");return EXIT_FAILURE;}
    ///options
    uint    collapseSquareW = 1,maxResSide=(uint) -1,unifyNearW = 0;
    ulong   maxOutSize=(ulong) -1,pixelDataLen,dataLen;
    char *ptr;
    const char* outFname=DFLT_OUTPATH;
    if (argc > 3){
        maxOutSize=strtoul(argv[3],&ptr,10);
        if (ptr==argv[3] ||  maxOutSize== ULONG_MAX){
            perror("maxOutSize  strtol errd");
            return EXIT_FAILURE;
        }
    }
    if (argc > 4){
        maxResSide =strtoul(argv[4],&ptr,10);
        if (ptr==argv[4] ||  maxResSide == UINT_MAX){
            perror("maxResSide strtol errd");
            return EXIT_FAILURE;
        }
    }
    if (argc > 5){
        collapseSquareW=strtoul(argv[5],&ptr,10);
        if (ptr==argv[5] ||  collapseSquareW == UINT_MAX){
            perror("collapseSquareW strtol errd");
            return EXIT_FAILURE;
        }
    }
    if (argc > 6)   outFname = argv[6];
    if (argc > 7){
        unifyNearW=strtoul(argv[7],&ptr,10);
        if (ptr==argv[7] ||  unifyNearW == UINT_MAX){
            perror("unifyNearW  strtol errd");
            return EXIT_FAILURE;
        }
    }
    data=malloc(sizeof(*data));
    if (!data){
        ERRPRINT("ppmData alloc for dense matrix failed\n");
        return EXIT_FAILURE;
    }
    data -> width  = ceil(N / collapseSquareW);
    data -> height = ceil(M / collapseSquareW);
    //checking max output size
    pixelDataLen = (N * M * PPM_CHAN_NUM ) / (collapseSquareW*collapseSquareW);
    if (pixelDataLen > maxOutSize){
        ERRPRINTS("THE GIVEN MATRIX IS EXCEEDING THE MAX SIZE OF %lu MB ",
            (pixelDataLen-maxOutSize) >>20 );
        collapseSquareW = sqrt(( (ulong) N * M * PPM_CHAN_NUM ) / maxOutSize);
        ERRPRINTS("SETTING collapseSquareW TO: %u\n",collapseSquareW);
        data -> width  = ceil(N / collapseSquareW);
        data -> height = ceil(M / collapseSquareW);
    }
    //checking max resolution longer side
    ulong matLongSide = MAX(data->width,data->height);
    if (matLongSide > maxResSide){
        ERRPRINTS("THE GIVEN MATRIX IS EXCEEDING MAX RESOLUTION SIDE OF %lu\n",
          matLongSide - maxResSide);
        collapseSquareW = matLongSide / maxResSide;
        data -> width  = ceil(N / collapseSquareW);
        data -> height = ceil(M / collapseSquareW);
        ERRPRINTS("SETTING collapseSquareW TO: %u\n",collapseSquareW);
    }   

    int headerLen = snprintf(data->header,PPM_HEADER_MAX_LEN,
        "P6\n%lu %lu\n255\n",data->width,data->height); 
    if (headerLen < 0){
        ERRPRINT("snprintf error");
        goto _free; 
    }
    if(__builtin_umull_overflow(data->height,data->width*PPM_CHAN_NUM,&pixelDataLen)){
        ERRPRINT("pixelDataLen overflow\n");
        goto _free;
    }
    if(__builtin_uaddl_overflow(pixelDataLen,headerLen,&dataLen)){
        ERRPRINT("dataLen overflow\n");
        goto _free;
    }

    printf("building ppm image in:%s, %lux%lu,\t with header:\n%s\n=>pixelDataLen: %luMB,"
      "collapseSquare:%u,unify:%u\n",outFname,data->width,data->height,
      data->header,dataLen>>20,collapseSquareW,unifyNearW);
    fflush(NULL);
    ///out file mmap for easy write
    outFd=open(outFname, O_RDWR | O_CREAT | O_EXCL | O_TRUNC, mode);
    if (errno==EEXIST)     outFd=open(outFname, O_RDWR | O_TRUNC, mode);
    if (outFd<0){
        perror("open outFd failed");
        goto _free;
    }
    if (ftruncate(outFd,dataLen)<0){
        perror("ftruncate err");
        goto _free;
    }
    map = mmap(NULL, dataLen, PROT_WRITE, MAP_SHARED, outFd, 0);    
    if (map == MAP_FAILED){
        perror("mmap failed... ");
        goto _free;
    }

    memcpy(map,data->header,headerLen);        //write header
    data->data=map+headerLen; //directly write converted matrix to outfile via mmap
    //memset(data->data,Z_PIXEL_COLOR,pixelDataLen); //ZERO=>IMPLICIT BY HOLE CREATED IN TRUNC
    printf("Alloc and prepare over, starting conversion\n");
    if (!strncmp(argv[2],MM_COO,strlen(MM_COO))){
        //denseMatrixToPPM(data,M,N,(double (*)[N]) denseMat,collapseSquareW,unifyNearW);
        if (sparseMatrixToPPM(data,sparseMat,collapseSquareW,unifyNearW))   goto _free;
    } 
    //TODO else PATTERN_DENSE
    
#ifdef TEST
    if (collapseSquareW!=1 || unifyNearW !=0)
        {ERRPRINT("TODO MAKE TEST CASE FOR THIS CONFIG");goto _free;}
    if (checkSparseMatrixToPPM(M,N, (uchar (*)[N]) data->data,
        sparseMat,collapseSquareW,unifyNearW))   goto _free;
    //TODO add support for step, unifyNearW, NOW CONSIDERED AS dflt... ->1a1 mat
    printf("DENSE MATCHING TEST PASSEED\n");
#endif
    out = EXIT_SUCCESS;
    
    _free:
    if(outFd)       close(outFd);
    if(data)        free(data);
    if(denseMat)    free(denseMat);
    if(map){
        if(munmap(map,dataLen) == -1){
            perror("Error un-mmapping the file");
        }
    //if(ftruncate(outFd,actualLen)<0){perror("ftruncate err ");goto _free;}//remove excess from mmapped
    }
    if (sparseMat)  freeSpmat(sparseMat);
    return out;
}
#endif
