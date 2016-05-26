// Threaded two-dimensional Discrete FFT transform
// Shivam Agarwal
// ECE4122 Project 2


#include <iostream>
#include <string>
#include <math.h>

#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.

Complex* ImageData;
int ImageWidth;
int ImageHeight;
int N = 1024;
Complex* Helper = new Complex[N*N];
Complex* W = new Complex[N/2];
int nThreads = 16;

// This is for MyBarrier
int P = nThreads + 1;
int count;
pthread_mutex_t countMutex;
bool* localCheck;
bool globalCheck;


unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// Call MyBarrier_Init once in main


void MyBarrier_Init()// you will likely need some parameters)
{
    count = P;
    
    pthread_mutex_init(&countMutex,0);
    localCheck = new bool[P];
    for(int i = 0; i < P; i++)
    {
        localCheck[i] = true;
    }
    globalCheck = true;
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier(int myId) // Again likely need parameters
{
    localCheck[myId] = !localCheck[myId];
    pthread_mutex_lock(&countMutex);
    int myCount = count;
    count--;
    pthread_mutex_unlock(&countMutex);
    if(myCount == 1)
    {
        count = P;
        globalCheck = localCheck[myId];
    }
    else
    {
        while(globalCheck != localCheck[myId]) { }
    }
}
                    
void Transform1D(Complex* h, int N)
{
  // Implementation of efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assuming even power of 2)
    Complex temp;
    Complex* W = new Complex[N/2];
    for(int n = 0; n < N/2; n++){
        W[n] = Complex(cos(2*M_PI*n/N), -sin(2*M_PI*n/N));
    }
    for(int np = 2; np <= N; np = np*2)
    {
        for(int i = 0; i < N; i = i + np)
        {
            for(int j = 0; j < np/2; j++)
            {
                int offset = np/2;
                temp = h[i + j];
                h[i + j] = h[i + j] + W[j*N/np] * h[i + j + offset];
                h[i + j + offset] = temp - W[j*N/np] * h[i + j + offset];
            }
        }
    }
}

void* Transform2DTHread(void* v)
{
    unsigned long myId = (unsigned long)v;
    
    int startRow = myId * N/nThreads;
    for(int n = 0; n < N/nThreads; n++)
    {
        int startElement = N*(startRow + n);
        Transform1D(ImageData + startElement, N);
    }
    
    
    MyBarrier(myId);
    
    return 0;
}

void Transform2D(const char* inputFN) 
{
    InputImage image(inputFN);
    ImageData = image.GetImageData();
    ImageWidth = image.GetWidth();
    ImageHeight = image.GetHeight();
    
    image.SaveImageData("MyInputImage.txt", ImageData, ImageWidth, ImageHeight);
    for(int n = 0; n < N; n++)
    {
        int startElement = n*N;
        Complex* temp = new Complex[N];
        for(int i = 0; i < N; i++)
        {
            temp[i] = ImageData[i + startElement];
        }
        for(int i = startElement; i < startElement + N; i++)
        {
            Helper[i] = temp[ReverseBits(i)];
        }
        delete [] temp;
    }
    
    for(int i = 0; i < N*N; i++)
    {
        ImageData[i] = Helper[i];
    }

    for(int n = 0; n < N/2; n++){
        W[n] = Complex(cos(2*M_PI*n/N), -sin(2*M_PI*n/N));
    }

    MyBarrier_Init();
    
    for (int i = 0; i < nThreads; ++i)
    {
        pthread_t pt;
        pthread_create(&pt, 0, Transform2DTHread, (void*)i);
    }
 
    MyBarrier(16);
    
    
    image.SaveImageData("MyAfter1D.txt", ImageData, ImageWidth, ImageHeight);
    
    int pos = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            Helper[pos] = ImageData[i + j*N];
            pos++;
        }
    }
    
    for(int i = 0; i < N*N; i++)
    {
        ImageData[i] = Helper[i];
    }

    for(int n = 0; n < N; n++)
    {
        int startElement = n*N;
        Complex* temp = new Complex[N];
        for(int i = 0; i < N; i++)
        {
            temp[i] = ImageData[i + startElement];
        }
        for(int i = startElement; i < startElement + N; i++)
        {
            Helper[i] = temp[ReverseBits(i)];
        }
        delete [] temp;
    }
    
    for(int i = 0; i < N*N; i++)
    {
        ImageData[i] = Helper[i];
    }
    
    MyBarrier_Init();
    
    for (int i = 0; i < nThreads; ++i)
    {
        pthread_t pt;
        pthread_create(&pt, 0, Transform2DTHread, (void*)i);
    }
    
    MyBarrier(16);

    pos = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            Helper[pos] = ImageData[i + j*N];
            pos++;
        }
    }
    
    for(int i = 0; i < N*N; i++)
    {
        ImageData[i] = Helper[i];
    }
    
    image.SaveImageData("MyAfter2D.txt", ImageData, ImageWidth, ImageHeight);
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  //if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  // MPI initialization here
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
