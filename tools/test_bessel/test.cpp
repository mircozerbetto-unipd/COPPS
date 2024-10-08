#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
extern "C"{
#include "bessel.h"
}



using namespace std;

int main (void)
{
	for (int i = 0; i < 100; i++)
		cout << (double)i*0.01 << "\t" << bessi(2,(double)i*0.01) << endl;
	return 0;
}
