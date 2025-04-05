#include "KSR_Task.h"
#include <iostream>
#include <locale>


int main()
{
	setlocale(LC_ALL, "Russian");
	double A{};
	double B{};
	double STEP{};
	double E_BORDER{};
	double E_ERROR{};
	int MAX_STEPS{};
	double START_POINT_FOR_U{};
	
	int option;
	
	std::cout << "1 - без ОЛП     2 - с ОЛП" << std::endl << std::endl;
	std::cout << "Выбирите вариант решения: ";
	std::cin >> option;
	std::cout << std::endl;

	double k;
	std::cout << "Введите параметр k: " << std::endl;
	std::cin >> k;

	double k_with_star;
	std::cout << "Введите параметр k*: " << std::endl;
	std::cin >> k_with_star;

	double c;
	std::cout << "Введите параметр c: " << std::endl;
	std::cin >> c;

	double m;
	std::cout << "Введите параметр m: " << std::endl;
	std::cin >> m;

	std::cout << "Введите левую границу интегрирования А: " << std::endl;
	std::cin >> A;

	std::cout << "Введите правую границу интегрирования B: " << std::endl;
	std::cin >> B;

	std::cout << "Введите шаг интегрирования STEP: " << std::endl;
	std::cin >> STEP;

	std::cout << "Введите E для оценки ЛП E_ERROR: " << std::endl;
	std::cin >> E_ERROR;

	std::cout << "Введите E-граничный E_BORDER: " << std::endl;
	std::cin >> E_BORDER;

	std::cout << "Введите максимальное число шагов MAX_STEPS: " << std::endl;
	std::cin >> MAX_STEPS;

	std::cout << "Введите начальное условие U(A): " << std::endl;
	std::cin >> START_POINT_FOR_U;


	
	KSR_Task Solution(A, B, STEP, E_ERROR, E_BORDER, MAX_STEPS, START_POINT_FOR_U);
	Solution.set_params(k, k_with_star, c, m);

	if (option == 1) {
		Solution.Solve_Without_Error_Control();
	}
	if (option == 2) {
		Solution.Solve_With_Error_Control();
	}

	//Solution.Write_To_File();
	Solution.PrintTable();
	Solution.PrintReference();

}