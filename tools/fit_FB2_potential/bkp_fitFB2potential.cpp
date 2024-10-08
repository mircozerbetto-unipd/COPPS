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
	if (argc < 8)
	{
		cout << endl << "Usage: " << endl;
		cout << "fitFB2potential xy-file xBins yBins xNmax yNmax min_x_angle_DEG min_y_angle_DEG" << endl << endl;
		return 1;
	}

	// 2. DEFINE SOME VARIABLES

	int xBins = atoi(argv[2]);
	int yBins = atoi(argv[3]);
	int xNmax = atoi(argv[4]);
	int yNmax = atoi(argv[5]);
	double minX = atof(argv[6]) * DEG2RAD;
	double minY = atof(argv[7]) * DEG2RAD;
	double maxX = minX + 2.0 * M_PI;
	double maxY = minY + 2.0 * M_PI;

	fstream iofile;

	// 3. PARSE INPUT FILE

	double xtmp, ytmp;
	dvector t1, t2;
	fstream infile;
	iofile.open(argv[1],ios::in);

	while (!iofile.eof())
	{
		iofile >> xtmp >> ytmp;
		t1.push_back(xtmp * DEG2RAD);
		t2.push_back(ytmp * DEG2RAD);
	}

	iofile.close();

	int nPoints = t1.size();

	// 4. MAKE 2D HISTOGRAM OF POTENTIAL
	
	cout << "* Making 2D histofgram of data..." << endl;
	bool count = false;

	double dx = (maxX - minX) / (double)(xBins);
	double dy = (maxY - minY) / (double)(yBins);

	dvector xEdges = dvector(xBins + 1, 0.0);
	dvector yEdges = dvector(yBins + 1, 0.0);

	for (int i = 0; i < xBins; i++) xEdges.at(i) = minX + (double)i * dx; xEdges.at(xBins) = maxX;
	for (int i = 0; i < yBins; i++) yEdges.at(i) = minY + (double)i * dy; yEdges.at(yBins) = maxY;

	dvector H = dvector(xBins * yBins, ZERO);

	for (int it = 0; it < nPoints; it++)
	{
		for (int ix = 0; ix < xBins; ix++)
		{
			count = t1.at(it) >= xEdges.at(ix) && t1.at(it) <= xEdges.at(ix + 1);
			if (count)
			{
				for (int iy = 0; iy < yBins; iy++)
				{
					count = t2.at(it) >= yEdges.at(iy) && t2.at(it) <= yEdges.at(iy + 1);
					if (count) H.at(ix * yBins + iy) += 1.0;
				}
			}
		}
	}
	
	double minU = 1.0e10;
	for (int i = 0; i < xBins; i++)
	{
		for (int j = 0; j < yBins; j++)
		{
			H.at(i * yBins + j) = -log(H.at(i * yBins + j));
			minU = H.at(i * yBins + j) < minU ? H.at(i * yBins + j) : minU;
		}
	}

	for (int i = 0; i < xBins * yBins; i++) H.at(i) -= minU;

	iofile.open("U_MD.dat",ios::out);
	for (int i = 0; i < xBins; i++)
	{
		for (int j = 0; j < yBins; j++)
			iofile << H.at(i * yBins + j) << "\t";
		iofile << endl;
	}
	iofile.close();
	cout << "* Created U_MD.dat file with the MD derived potential of mean force" << endl;

	// 5. FIND WEIGHTS OF THE EXPANSION

	cout << "* Projecting the potential on the basis..." << endl;
	int dim1 = 2 * xNmax + 1;
	int dim2 = 2 * yNmax + 1;
	double a1, a2, angle, H12;
	double norm = dx * dy / (4.0 * M_PI * M_PI);
	cdouble czero = cdouble(0.0, 0.0);
	cvector weights = cvector(dim1 * dim2, czero);

	for (int ix = 0; ix < xBins; ix++)
	{
	        a1 = 0.5 * (xEdges[ix] + xEdges[ix+1]);
	        for (int iy = 0; iy < yBins; iy++)
		{
	                a2 = 0.5 * (yEdges[iy] + yEdges[iy+1]);
	                H12 = norm * H.at(ix * yBins + iy);
	                for (int n1 = -xNmax; n1 <= xNmax; n1++)
			{
                        	for (int n2 = -yNmax; n2 <= yNmax; n2++)
				{
					angle = (double)n1 * a1 + (double)n2 * a2;
					weights[(n1 + xNmax) * dim2 + (n2 + yNmax)] -= H12 * cdouble(cos(angle), sin(angle));
				}
			}
		}
	        cout  << "* " << ix+1 << " / " << xBins << endl;
	}
	cout << "* ...projections done" << endl;

	// 6. CALCULATE THE POMF
	
	cout << "* Calculateing the fitted potential on initial grid..." << endl;
	double m;
	cvector U = cvector(xBins * yBins, czero);
	minU = 1.0e10;
	for (int ix = 0; ix < xBins; ix++)
	{
		a1 = 0.5 * (xEdges[ix] + xEdges[ix+1]);
		for (int iy = 0; iy < yBins; iy++)
		{
			a2 = 0.5 * (yEdges[iy] + yEdges[iy+1]);
			U.at(ix * yBins + iy) = 0.0;
			for (int n1 = -xNmax; n1 <= 0; n1++)
			{
				int n2UpLim = n1 < 0 ? yNmax : 0;
				for (int n2 = -yNmax; n2 <= n2UpLim; n2++)
				{
					m = n1 == 0 && n2 == 0 ? 0.5 : 1.0;
					angle = (double)n1 * a1 + (double)n2 * a2;
					U.at(ix * yBins + iy) -= m * ( weights[(n1 + xNmax) * dim2 + (n2 + yNmax)] * cdouble(cos(angle), -sin(angle)) + weights[(-n1+xNmax)*dim2+(-n2+yNmax)] * cdouble(cos(angle),sin(angle)));
				}

			}
			minU = U.at(ix * yBins + iy).real() < minU ? U.at(ix * yBins + iy).real() : minU;
		}
	        cout << "* " << ix+1 << " / " << xBins << endl;
	}

	// Reset the minimum of U to 0
	for (int i = 0; i < xBins * yBins; i++) U.at(i) -= minU;
	weights.at(xNmax * dim2 + yNmax) -= minU;

	cout << "* ...potential of mean force ready." << endl;

	iofile.open("POMF.dat",ios::out);
	for (int i = 0; i < xBins; i++)
	{
		for (int j = 0; j < yBins; j++)
			iofile << U.at(i * yBins + j).real() << "\t";
		iofile << endl;
	}
	iofile.close();
	cout << "* Created POMF.dat file with the fitted potential of mean force" << endl;

	// 7. MEASURE CHI-SQUARE BETWENN SURFACES

	double chisquare = 0.0;
	for (int i = 0; i < xBins * yBins; i++)
                chisquare += (U.at(i).real() - H.at(i)) * (U.at(i).real() - H.at(i));
	chisquare /= (double)(xBins * yBins);
	cout << "* Reduced Chi-square = " << chisquare << endl;

	// 8. OUTPUT FOR C++OPPS 2.0

	iofile.open("fb2_pomf_copps.dat",ios::out);
	iofile << "potential_coefficients_t2:" << xNmax << ":" << yNmax << endl;
	for (int n1 = -xNmax; n1 <= xNmax; n1++)
	{
		for (int n2 = -yNmax; n2 <= yNmax; n2++)
			iofile << weights[(n1 + xNmax) * dim2 + (n2 + yNmax)].real() << " " << weights[(n1 + xNmax) * dim2 + (n2 + yNmax)].imag() << " ";
		iofile << endl;
	}
	iofile.close();
	cout << "* Created fb2_pomf_copps.dat file with the input to include in _copps.input (v >= 2.0) file to use this POMF" << endl;
	cout << "* IMPORTANT: POMF coefficients are given in kT units at the temperature at which the molecular dynamics trajectory was calculated"<<endl;

	// 9. OUTPUT FILE FOR VISUALIZATION WITH OCTAVE
	
	iofile.open("plot_pomf.m",ios::out);

	iofile << "x = [";
	for (int i = 0; i < xBins; i++) iofile << 0.5 * (xEdges[i] + xEdges[i+1]) * RAD2DEG << " ";
	iofile << "];" << endl;
		iofile << "y = [";
	for (int i = 0; i < yBins; i++) iofile << 0.5 * (yEdges[i] + yEdges[i+1]) * RAD2DEG << " ";
	iofile << "];" << endl;
	iofile << "[X,Y] = meshgrid(x,y);" << endl;

	iofile << "figure(1);" << endl;
	iofile << "Z1 = load('./U_MD.dat');" << endl;
	iofile << "s1 = surf(X,Y,Z1); set(s1, 'edgecolor', 'none'); view(0,90); \%colorbar;" << endl;
	iofile << "xlabel('t1 / deg'); ylabel('t2 / deg');" << endl;

	iofile << "figure(2);" << endl;
	iofile << "Z2 = load('./POMF.dat');" << endl;
	iofile << "s2 = surf(X,Y,Z2); set(s2, 'edgecolor', 'none'); view(0,90); \%colorbar;" << endl;
	iofile << "xlabel('t1 / deg'); ylabel('t2 / deg');" << endl;

	return 0;
}
