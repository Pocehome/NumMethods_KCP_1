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
	
	std::cout << "1 - ��� ���     2 - � ���" << std::endl << std::endl;
	std::cout << "�������� ������� �������: ";
	std::cin >> option;
	std::cout << std::endl;

	double k;
	std::cout << "������� �������� k: " << std::endl;
	std::cin >> k;

	double k_with_star;
	std::cout << "������� �������� k*: " << std::endl;
	std::cin >> k_with_star;

	double c;
	std::cout << "������� �������� c: " << std::endl;
	std::cin >> c;

	double m;
	std::cout << "������� �������� m: " << std::endl;
	std::cin >> m;

	std::cout << "������� ����� ������� �������������� �: " << std::endl;
	std::cin >> A;

	std::cout << "������� ������ ������� �������������� B: " << std::endl;
	std::cin >> B;

	std::cout << "������� ��� �������������� STEP: " << std::endl;
	std::cin >> STEP;

	std::cout << "������� E ��� ������ �� E_ERROR: " << std::endl;
	std::cin >> E_ERROR;

	std::cout << "������� E-��������� E_BORDER: " << std::endl;
	std::cin >> E_BORDER;

	std::cout << "������� ������������ ����� ����� MAX_STEPS: " << std::endl;
	std::cin >> MAX_STEPS;

	std::cout << "������� ��������� ������� U(A): " << std::endl;
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