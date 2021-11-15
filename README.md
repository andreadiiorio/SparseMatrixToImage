CONVERT SPARSE MATRIX MARKET COO INTO A PPM image
=================================================

Convert a sparse matrix in matrix market coordinate format into a PPM image
(easy to compress using cjpeg tool)

the output file will be mmapped, the 0 filling is exploited for the ZERO elements, 
the NON ZERO elements will be written with a with white (RGB pixels 255,255,255) dot.

Along with the sparse matrix to image logic it's shipped also some utils and macros 
for the reading and parsing of the sparse matrix (from my works SPGEMV, SPGEMM FOR AMG) //TODO LINK HERE

#USAGE
COMPILE JUST TYPING MAKE

USAGE:
./sparseMatrixToImage.o <inputMatrixMarket_Coo_file, MM_COO> [MAX_OUT_SIZE, collapseSquareW, outFilePath]

where: 
-collapseSquareW -> collapse a 2D square of the given width into a single dot in the output img
-MAX_OUT_SIZE    -> Maximum size for the output file, if the matrix dimension * 3 (colors channels num)
                    exceed this threshold, the collapseSquareW will be enlarged to
                    sqrt( (N*M*3)/MAX_OUT_SIZE )

