int counter=-1, nSteps_;

double xPosi_[7000];
double yPosi_[7000];
double zPosi_[7000];
double intensity_[4][7000] ={{0.0}};
int MaxEddy_ = 7000;
double fsigma_[7000] = {0.0};


double randomArray[7000] = {0}, randomArray1[7000] = {0}, randomArray2[7000] =  {0}, randomArray3[7000] =  {0} , randomArray4[7000] =  {0};


double RST11 =0, RST22 =0 , RST33 = 0, RST21 = 0, RST31 = 0, RST32 = 0;
double a11[500][500] = {{0.0}}, a22[500][500] = {{0.0}}, a33[500][500] = {{0.0}}, a21[500][500] = {{0.0}}, a31[500][500] = {{0.0}}, a32[500][500] = {{0.0}};


int indexI[25000] = {0};
int indexJ[25000] = {0};

double Ustored[300] = {0.0};
double Vstored[300] = {0.0};
double Wstored[300] = {0.0};
double storedRad[300] = {0.0};

double upert[500][500] = {{0.0}}, vpert[500][500] = {{0.0}}, wpert[500][500] = {{0.0}};


double sigma_[500][500]={{0.0}};

double sigmaMax_= -1 , sigmaMin_ =100;



double a12=0.0;
double a13=0.0;
double a23=0.0;
