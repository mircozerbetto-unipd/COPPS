#include "redwig.h"

using namespace std;

int main(void)
{
	int j=0,m,k;
	long double beta;
	wigner w(30);
	while (j >=0)
	{
		cout << "Input j,m,k: ";
		cin >> j >> m >> k;
		cout << "Input beta / deg: ";
		cin >> beta;
		beta *= M_PI/180.0;
		cout << "d = " << w.getReducedWignerMatrix(j,m,k,beta) << endl;
	}
	return (0);
}
