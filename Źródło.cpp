#include <iostream>
#include <ctime>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <random>

using namespace std;


double function(double x[], int n, int eq);
double sphere_function(double x[], int n);
double levy_function(double x[], int n);
double himmelblau_function(double x[]);
double ackley_function(double x[]);
double rastrigin_function(double x[]);
double rosenbrock_function(double x[]);
double schwefel_function(double x[]);
double michalewicz_function(double x[]);
random_device R;

int f_counter = 0;	//licznik wywo³añ funkcji celu da tej próby
int g_counter = 0;	//licznik wywo³añ funkcji celu dla wszystkich prób
double outcome = 0; // suma wyników najlepszych funkcji celu po osi¹gniêciu warunku stopu
int stop_cnd = 5000; // maksymalna iloœc wywo³ania funkcji
double y_eps;
bool reset = false;


class ant {
public:
	double* x;
	double y;
};
class colony {
public:
	int eq = 1;
	int n = 10;
	int m = 10;
	int ant_amount = 10;
	int dim = 1;
	double min = -10, max = 10;
	ant* ants = new ant[ant_amount];
	int* pheromone = new int[ant_amount]();
	void create_ants(double min, double max)	//losuje mrówki i zapisuje si³e feromonów w tablicy pheromone od najni¿szego wyniku
	{
		for (int i = 0; i < ant_amount; ++i)
		{
			ant an;
			an.x = new double[dim];
			for (int i = 0; i < dim; ++i)
			{
				double f = 1.0 * R() / R.max();
				an.x[i] = (min)+f * (max - min);
			}
			an.y = function(an.x, dim, eq);
			ants[i] = an;

			int j = 0;
			int num = pheromone[j];
			while (ants[num].y < ants[i].y)
			{
				++j;
				num = pheromone[j];
				if (j >= i)
					break;
			};
			for (int k = i; k > j; --k)
			{
				pheromone[k] = pheromone[k - 1];
			}
			pheromone[j] = i;
		}
	}
	void move(double a)
	{
		for (int i = ant_amount - 1; 0 < i; --i)
		{
			int num = pheromone[i];
			double distance[2];
			int num2 = pheromone[i - 1];
			double suma = 0;
			for (int i = 0; i < dim; ++i)
			{
				suma += pow(ants[num].x[i] + ants[num2].x[i], 2);
			}
			distance[0] = sqrt(suma);	//zapisanie pierwszej odleg³oœci 
			distance[1] = i - 1;	//zapisanie punktu do ktorego mierzylismy odleglosc( zapisane w id feromonu nie id mrówki)
			for (int j = i - 1; 0 < j; --j)
			{
				suma = 0;
				num2 = pheromone[j - 1];
				//sprawdzamy która mrówka z silniejszym feromonem jest najbli¿ej

				for (int i = 0; i < dim; ++i)
				{
					suma += pow(ants[num].x[i] + ants[num2].x[i], 2);
				}
				double check = sqrt(suma);
				if (check < distance[0] || distance[0] == 0)
				{
					distance[0] = check;
					distance[1] = j - 1;
				}
			}
			//zaczynamy przesuwanie
			double f = 1.0 * R() / R.max();
			double mv = (0)+f * (a/10);
			int pheromone_id = distance[1];
			num2 = pheromone[pheromone_id];//mrówka do której mrówka któr¹ przesuwamy ma najbli¿ej i która ma silniejszy œlad feromonowy
			double dis;
			int j = 0;
			int k = 0;
			double tmp;
			for (j; j < dim; ++j)	//przesuniêcie po pierwszej wspó³¿êdnej która mo¿e byæ zmieniona
			{
				if (ants[num].x[j] == ants[num2].x[j])
					continue;
				else
				{
					if (ants[num].x[j] > ants[num2].x[j])
						mv = -mv;
					dis = abs(ants[num].x[j] - ants[num2].x[j]);	// zapisany dystans na danej osi miêdzy mrówk¹ która idzie a punktem do którego siê kieruje
					if ((ants[num].x[j] + mv)< max && (ants[num].x[j] + mv)> min)
						ants[num].x[j] = ants[num].x[j] + mv;
					else
						continue;
					k = j;
					++j;
					break;
				}
			}
			double wsp = abs(mv / dis);
			for (j;j < dim;++j)
			{
				if (ants[num].x[j] == ants[num2].x[j])
					continue;
				else
				{
					tmp = ants[num].x[j] + (ants[num2].x[j] - ants[num].x[j]) * wsp;
					if(tmp < max && tmp > min)
					ants[num].x[j] = ants[num].x[j] + (ants[num2].x[j] - ants[num].x[j]) * wsp;
					else
						ants[num].x[j] = ants[num].x[j] - mv;
				}
			}
			num = pheromone[i];
			ants[num].y = function(ants[num].x, dim, eq);
			if (reset)
				return;
		}
		//mrówki o feromonach 1-9 zosta³y przesuniête w kierunku najbli¿szej lepszej mrówki
		//teraz wybieramy kierunek i przesuwamy najlepsz¹ mrówke;

		bool* direction = new bool[dim];
		int num = pheromone[0];
		ant tmp;
		tmp.x = new double[dim];
		//double res;
		for (int i = 0; i < dim; ++i)
			tmp.x[i] = ants[num].x[i];
		for (int i = 0; i < dim; ++i)
		{
			tmp.x[i] = tmp.x[i] * 1.1;
			if (reset)
				return;
			if (function((tmp.x), dim, eq) < (function((ants[num].x), dim, eq)))
				direction[i] = false;
			else
				direction[i] = true;
			if (reset)
				return;
			if (ants[num].x < 0)
				direction[i] = !direction[i];
			tmp.x[i] = ants[num].x[i];
			double f = 1.0 * R() / R.max();
			double mv = (0)+f * (a/10);
			if (!direction[i])
				mv = -mv;//zmiana kireunku przesuniecia je¿eli direction jest false

			if ((ants[num].x[i] + mv) < max && (ants[num].x[i] + mv) > min)
			{
				tmp.x[i] = tmp.x[i] + mv;
				if (function((tmp.x), dim, eq) < (function((ants[num].x), dim, eq)))
					ants[num].x[i] = tmp.x[i];
				else
					tmp.x[i] = ants[num].x[i];
				if (reset)
					return;
			}
			else
				tmp.x[i] = ants[num].x[i];
		}
		ants[num].y = (function((ants[num].x), dim, eq));
	}
	void quick(int* pheromone, int start, int end)
	{
		double number = ants[pheromone[start]].y;
		int mniejsza = start, wieksze = start + 1;
		for (int i = start + 1; i <= end; ++i)
		{
			if (ants[pheromone[i]].y < number)
			{
				swap(pheromone[i], pheromone[wieksze]);
				++mniejsza;
				++wieksze;
			}
		}
		swap(pheromone[start], pheromone[mniejsza]);
		if(start!=mniejsza)
		quick(pheromone, start, mniejsza - 1);
		if(mniejsza!=end)
		quick(pheromone, mniejsza + 1, end);
	}

};

double function(double x[], int n, int eq)
{
	++f_counter;
	double result;
	switch (eq)
	{
	case 1:
		result = sphere_function(x, n);
		break;
	case 2:
		result = michalewicz_function(x);
		break;
	case 3:
		result = schwefel_function(x);
		break;
	case 4:
		result = ackley_function(x);
		break;
	case 5:
		result = rastrigin_function(x);
		break;
	case 6:
		result = rosenbrock_function(x);
		break;
	case 7:
		result = levy_function(x, n);
		break;
	}
	if (f_counter >= stop_cnd || result <= y_eps)
		reset = true;
	return result;
}
double sphere_function(double x[], int n)
{
	double result = 0;
	for (int i = 0; i < n;++i)
	{
		result += pow(x[i], 2);
	}
	return result;
}
double levy_function(double x[], int n)
{
	double result = 0;
	for (int i = 0; i < n;++i)
	{
		result += pow((1 + (x[i] - 1) / 4) - 1, 2) * (1 + 10 * pow(sin(M_PI * (1 + (x[i] - 1) / 4) + 1),2));
	}
	result += pow(sin(M_PI * (1 + (x[0] - 1) / 4)),2) + pow((1 + (x[n - 1] - 1) / 4) - 1, 2) * (1 + 10 * pow(sin(M_PI * (1 + (x[n - 1] - 1) / 4) + 1), 2));
	return result;
}
double michalewicz_function(double x[])
{
	double result;
	result = -sin(x[0]) * pow(sin(pow(x[0], 2) / M_PI), 20) - sin(x[1]) * pow(sin(2 * pow(x[1], 2) / M_PI), 20);
	return result;
}
double schwefel_function(double x[])
{
	double result;
	result = (-x[0]) * sin(sqrt(abs(x[0]))) - (x[1] * sin(sqrt(abs(x[1]))));
	return result;
}
double ackley_function(double x[])
{
	double result;
	result = -20 * exp(-0.2 * sqrt(0.5 * (pow(x[0], 2) + pow(x[1], 2)))) - exp(0.5 * (cos(2 * M_PI * x[0]) + cos(2 * M_PI * x[1]))) + M_E + 20;
	return result;
}
double rastrigin_function(double x[])
{
	double result;
	result = pow(x[0], 2) + pow(x[1], 2) - cos(18 * x[0]) - cos(18 * x[1]) + 2;
	return result;
}
double rosenbrock_function(double x[])
{
	double result;
	result = 100 * (pow(x[1] - pow(x[0], 2), 2)) + pow((1 - x[0]), 2);
	return result;
}
int main()
{
	int tests_number = 100;
	for (int i = 0; i < tests_number; ++i)
	{
		reset = false;
		int eq = 7; // wybór dzia³ania
		int n = 10;
		int m = 30;
		int ants = 10;
		double min = -5;
		double max = 15;
		double a = (max - min) / 5;
		stop_cnd = 5000;
		y_eps = pow(10, -3);
		colony test;
		test.dim = 5;
		test.eq = eq;
		test.ant_amount = ants;
		test.ants[ants];
		test.create_ants(min, max);
		test.min = min;
		test.max = max;

		for (int i = 0; i < n; ++i)
		{
			for (int i = 0; i < m;++i)
			{
				test.move(a);
				test.quick(test.pheromone, 0, ants - 1);
				if (reset)
					break;
			}
			if (reset)
				break;
			a = a / 10;
		}


		double best = test.ants[0].y;
		for (int i = 1; i < test.ant_amount; ++i)
		{
			if (test.ants[i].y < best)
				best = test.ants[i].y;
		}
		outcome += best;
		g_counter += f_counter;
		f_counter = 0;
	}

	cout << endl << endl << "*********+++++++++++++**********" << endl;
	cout << "ilosc wywolan funkcji = " << g_counter/tests_number<< endl;
	cout << "srednia odnalezina wartosc = " << outcome/tests_number<< endl;
	return 0;
}