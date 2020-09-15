namespace imajuscule::bmp {

/*
int main ()
{
  int height = 361;
  int width = 867;
  unsigned char image[height][width][BYTES_PER_PIXEL];
  char* imageFileName = (char*) "bitmapImage.bmp";
  
  int i, j;
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      image[i][j][2] = (unsigned char) ( i * 255 / height );             ///red
      image[i][j][1] = (unsigned char) ( j * 255 / width );              ///green
      image[i][j][0] = (unsigned char) ( (i+j) * 255 / (height+width) ); ///blue
    }
  }
  
  generateBitmapImage((unsigned char*) image, height, width, imageFileName);
  printf("Image generated!!");
}
*/

void generateBitmapImage (unsigned char const* image, int height, int width, char const* imageFileName);

} // NS
