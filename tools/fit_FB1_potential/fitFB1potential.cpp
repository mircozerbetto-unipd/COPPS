#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

#define ZERO 1.0e-10
#define DEG2RAD M_PI / 180.0
#define RAD2DEG 180.0 / M_PI

typedef vector <double> dvector;

typedef complex <double> cdouble;
typedef vector <cdouble> cvector;

int main (int argc, char* argv[])
{
	// 1. CHECK INPUT

	if (argc < 6)
	{
		cout << endl << "Usage: " << endl;
		cout << "fitFB1potential 1-column-file nBins nMax min_angle_DEG U_flag" << endl << endl;
		return 1;
	}

	// 2. DEFINE SOME VARIABLES

	int nBins = atoi(argv[2]);
	int nMax = atoi(argv[3]);
	double minX = atof(argv[4]) * DEG2RAD;
	int haveU = atoi(argv[5]);
	double maxX = minX + 2.0 * M_PI;

	fstream iofile;

	double dx = (maxX - minX) / (double)(nBins);
	dvector xEdges = dvector(nBins + 1, 0.0);
	for (int i = 0; i < nBins; i++) xEdges.at(i) = minX + (double)i * dx; xEdges.at(nBins) = maxX;
	dvector H = dvector(nBins, ZERO);
	double minU = 1.0e10;

	if (!haveU)
	{

		// 3. PARSE INPUT FILE

		double xtmp;
		dvector t1;
		iofile.open(argv[1],ios::in);

		while (!iofile.eof())
		{
			iofile >> xtmp;
			t1.push_back(xtmp * DEG2RAD);
		}

		iofile.close();

		int nPoints = t1.size();

		// 4. MAKE 1D HISTOGRAM OF POTENTIAL
		
		cout << "* Making 1D histofgram of data..." << endl;
		bool count = false;

		for (int it = 0; it < nPoints; it++)
		{
			for (int ix = 0; ix < nBins; ix++)
			{
				count = t1.at(it) >= xEdges.at(ix) && t1.at(it) <= xEdges.at(ix + 1);
				if (count) H.at(ix) += 1.0;
			}
		}
		
		for (int i = 0; i < nBins; i++)
		{
			H.at(i) = -log(H.at(i));
			minU = H.at(i) < minU ? H.at(i) : minU;
		}
	}
	else // U must already be in kT units
	{
		iofile.open(argv[1],ios::in);
		for (int ix = 0; ix < nBins; ix++)
		{
			iofile >> H.at(ix);
			minU = H.at(ix) < minU ? H.at(ix) : minU;
		}
		iofile.close();
	}

	for (int i = 0; i < nBins; i++) H.at(i) -= minU;

	iofile.open("U_MD.dat",ios::out);
	for (int i = 0; i < nBins; i++)
	{
		iofile << H.at(i) << endl;
	}
	iofile.close();
	cout << "* Created U_MD.dat file with the MD derived potential of mean force" << endl;

	// 5. FIND WEIGHTS OF THE EXPANSION

	cout << "* Projecting the potential on the basis..." << endl;
	int dim1 = nMax + 1;
	double a1, H1;
	double norm = dx / (2.0 * M_PI);
	cdouble czero = cdouble(0.0, 0.0);
	cvector weights = cvector(dim1, czero);

	for (int ix = 0; ix < nBins; ix++)
	{
		a1 = 0.5 * (xEdges[ix] + xEdges[ix+1]);
		H1 = norm * H.at(ix);
		for (int n1 = 0; n1 <= nMax; n1++)
			weights[n1] -= H1 * cdouble(cos((double)n1 * a1), sin((double)n1 * a1));
	}
	cout << "* ...projections done" << endl;

	// 6. CALCULATE THE POMF
	
	cout << "* Calculateing the fitted potential on initial grid..." << endl;
	double m;
	dvector U = dvector(nBins, ZERO);
	minU = 1.0e10;
	for (int ix = 0; ix < nBins; ix++)
	{
		a1 = 0.5 * (xEdges[ix] + xEdges[ix+1]);
		U.at(ix) = 0.0;
		for (int n1 = 0; n1 <= nMax; n1++)
		{
			m = n1 == 0 ? 1.0 : 2.0;
			U.at(ix) -= m * ( weights.at(n1).real() * cos((double)n1 * a1) + weights.at(n1).imag() * sin((double)n1 * a1) );
		}
		minU = U.at(ix) < minU ? U.at(ix) : minU;
	}

	// Reset the minimum of U to 0
	for (int i = 0; i < nBins; i++) U.at(i) -= minU;
	weights.at(nMax) -= minU;

	cout << "* ...potential of mean force ready." << endl;

	iofile.open("POMF.dat",ios::out);
	for (int i = 0; i < nBins; i++)
		iofile << U.at(i) << endl;
	iofile.close();
	cout << "* Created POMF.dat file with the fitted potential of mean force" << endl;

	// 7. MEASURE CHI-SQUARE BETWENN SURFACES

	double chisquare = 0.0;
	for (int i = 0; i < nBins; i++)
                chisquare += (U.at(i) - H.at(i)) * (U.at(i) - H.at(i));
	chisquare /= (double)(nBins);
	cout << "* Reduced Chi-square = " << chisquare << endl;

	// 8. OUTPUT FOR C++OPPS 2.0

	iofile.open("fb1_pomf_copps.dat",ios::out);
	for (int n1 = 0; n1 <= nMax; n1++)
		iofile << "potential_coefficient_t1:" << weights.at(n1).real() << ":" << weights.at(n1).imag() << ":" << "kT" << endl;
	iofile.close();
	cout << "* Created fb1_pomf_copps.dat file with the input to include in _copps.input (v >= 2.0) file to use this POMF" << endl;
	cout << "* IMPORTANT: POMF coefficients are given in kT units at the temperature at which the molecular dynamics trajectory was calculated"<<endl;

	// 9. OUTPUT FILE FOR VISUALIZATION WITH OCTAVE
	
	iofile.open("plot_pomf.m",ios::out);

	iofile << "x = [";
	for (int i = 0; i < nBins; i++) iofile << 0.5 * (xEdges[i] + xEdges[i+1]) * RAD2DEG << " ";
	iofile << "];" << endl;

	iofile << "u1 = load('./U_MD.dat');" << endl;
	iofile << "u2 = load('./POMF.dat');" << endl;
	iofile << "plot(x,u1,x,u2);" << endl;
	iofile << "xlabel('t1 / deg');" << endl;

	iofile.close();

	return 0;
}
