#include "TestTask.h"
#include "TrueSolution.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <queue>

#include <iomanip>
#include <iostream>


//��������������� �������
bool TestTask::x_in_border(const double& B, const double& X, const double& BORDER) { // ������� ��������� , ��������� �� ������� X � ����������� ������ �������
	return X >= B - BORDER && X <= B;
}



//�������� �������
TestTask::TestTask(const double& _A, const double& _B, const double& _STEP, const double& _E_ERROR, const double& _E_BORDER, const int& _MAX_STEPS, const double& _START_POINT) {
	parametrs.A = _A;
	parametrs.B = _B;
	parametrs.STEP = _STEP;
	parametrs.E_ERROR = _E_ERROR;
	parametrs.E_BORDER = _E_BORDER;
	parametrs.MAX_STEPS = _MAX_STEPS;
	parametrs.START_POINT_FOR_U = _START_POINT;
}

TestTask::~TestTask(){
	parametrs.A = 0.0;
	parametrs.B = 0.0;
	parametrs.STEP = 0.0;
	parametrs.E_ERROR = 0.0;
	parametrs.E_BORDER = 0.0;
	parametrs.MAX_STEPS = 0;
	parametrs.START_POINT_FOR_U = 0.0;

	TABLE_INFORMATION.clear();
	DISTANCE_Ui_Vi.clear();
	STEPS_and_Xs.clear();
	ERRORS_LIST.clear();
}


double TestTask::f(const double& X, const double& V)
{
	return (-1)*(V/2);
}

double TestTask::find_K1(const double& X, const double& V)
{
	return f(X,V);
}

double TestTask::find_K2(const double& X, const double& V, const double& STEP, const double& K1)
{
	return f(X+STEP/2.0,V+(STEP/2.0)*K1);
}

double TestTask::find_K3(const double& X, const double& V, const double& STEP, const double& K2)
{
	return f(X + STEP / 2.0, V + (STEP / 2.0) * K2);
}

double TestTask::find_K4(const double& X, const double& V, const double& STEP, const double& K3)
{
	return f(X + STEP, V + STEP * K3);
}

void TestTask::make_Step(double& X, double& V, const double& STEP)
{
	double K1 = find_K1(X, V);
	double K2 = find_K2(X, V, STEP, K1);
	double K3 = find_K3(X, V, STEP, K2);
	double K4 = find_K4(X, V, STEP, K3);
	X += STEP;
	V += (STEP / 6.0) * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
}

int TestTask::control_Error(double& X,double& V, double& X_EXTRA, double& V_EXTRA, double& OLD_X, double& OLD_V, double& S, double& CURRENT_DOUBLING, double& CURRENT_REDUCTION)
{
	if (abs(S) < parametrs.E_ERROR / pow(2, 5)) { // ������� ��� ��������
		parametrs.STEP *= 2;
		++reference.STEP_DOUBLING_COUNT;
		++CURRENT_DOUBLING;

		return 0;
	}
	else if (abs(S) > parametrs.E_ERROR) {        //���� ������� ��� ������� ���� �����������, �� ������������ �����,
		bool FLAG_TO_EXIT = false;
		double h = parametrs.STEP;  				  // ��������� ������ � ���������� �����, ���� ��� ���� ����������� ����� 
		while (!FLAG_TO_EXIT) {	
			                                 // ��������� �������, �� ����� ������������ � ����� ���, � ��� ���� ����������� �� ����� ����������
			h /= 2;
			++CURRENT_REDUCTION;
			++reference.STEP_REDUCTION_COUNT;

			V = OLD_V;
			X = OLD_X;
			V_EXTRA = V;
			X_EXTRA = X;

			for (int i = 0; i < 2; ++i) {
				make_Step(X_EXTRA, V_EXTRA, h / 2.0);
			}
			make_Step(X, V, h);

			S = (V_EXTRA - V) / (pow(2, 4) - 1);
			if (abs(S) <= parametrs.E_ERROR) {
				FLAG_TO_EXIT = true;
				parametrs.STEP = h;
			}
		}

		return 1;
	}
}

void TestTask::find_Max_Step()
{
	reference.MAX_STEP_WITH_X = *std::max_element(STEPS_and_Xs.begin(), STEPS_and_Xs.end(),
		[](const std::pair<double, double>& a, const std::pair<double, double>& b) {
			return a.first < b.first;
		});
}

void TestTask::find_Min_Step()
{
	reference.MIN_STEP_WITH_X = *std::min_element(STEPS_and_Xs.begin(), STEPS_and_Xs.end(),
		[](const std::pair<double, double>& a, const std::pair<double, double>& b) {
			return a.first < b.first;
		});
}

void TestTask::find_Max_Error() {
	reference.MAX_ERROR = *std::max_element(ERRORS_LIST.begin(), ERRORS_LIST.end());
}

void TestTask::find_Min_Error() {
	reference.MIN_ERROR = *std::min_element(ERRORS_LIST.begin(), ERRORS_LIST.end());
}

void TestTask::find_Last_X() {
	reference.LAST_X = TABLE_INFORMATION.back()[1];
}

void TestTask::find_Last_V() {
	reference.LAST_V = TABLE_INFORMATION.back()[2];
}

void TestTask::Solve_Without_Error_Control()
{
	double X = parametrs.A;
	double U = find_True_Solution(X,parametrs.START_POINT_FOR_U); //������� �������� ������� �� ������� ����
	double V = U;
	
	std::vector<double> TABLE_ROW1 = {0.0, X, V, parametrs.STEP, U, abs(U-V)}; // ����� � ����� - ������ �������� ������� � ���� { i X_i V_i STEP_i U_i |U_i - V_i| }
	TABLE_INFORMATION.emplace_back(TABLE_ROW1);

	double OLD_X = X; //� ��������� OLD ������ �������� � ���������� ����
	double OLD_V = V;
	double OLD_U = U;

	bool FLAG_TO_EXIT = false;
	for (int i = 1; i <= parametrs.MAX_STEPS; ++i)
	{
		make_Step(X, V, parametrs.STEP);
		U = find_True_Solution(X, parametrs.START_POINT_FOR_U);

		if (X > parametrs.B + parametrs.E_BORDER) {           
													// ���� X ����� �� ������ �������, ������������ �� ��� ����� � ������ ���,                                        
			X = OLD_X;								//  ������ ������� ������ ������� � ��������� ������ X � ����������� ��������������
			V = OLD_V;
			U = OLD_U;
			parametrs.STEP = parametrs.B - OLD_X;
			make_Step(X, V, parametrs.STEP);
			U = find_True_Solution(X, parametrs.START_POINT_FOR_U);


			FLAG_TO_EXIT = true;                        //���� X ������ � ������ ��������, ����������� �������������� 
			
		}
		if (x_in_border(parametrs.B, X, parametrs.E_BORDER)) {
			FLAG_TO_EXIT = true;
		}

		DISTANCE_Ui_Vi.emplace_back(abs(U - V));
		std::vector<double> TABLE_ROW = { static_cast<double>(i), X, V, parametrs.STEP, U, abs(U - V)};
		TABLE_INFORMATION.emplace_back(TABLE_ROW);
		++reference.ITERATIONS_COUNT;

		if (FLAG_TO_EXIT)
			break;
		else {
			OLD_X = X;
			OLD_V = V;
			OLD_U = U;
		}

		
	}
	find_Max_Distance_Ui_Vi();
}

void TestTask::Solve_With_Error_Control()
{
	double X = parametrs.A;
	double U = find_True_Solution(X, parametrs.START_POINT_FOR_U);
	double V = U;

	std::vector<double> TABLE_ROW1 = { 0.0, X, V, V, 0.0, 0.0, parametrs.STEP, 0.0, 0.0, U, abs(U - V) }; // ����� � ����� - ������ �������� ������� � ���� 
	TABLE_INFORMATION.emplace_back(TABLE_ROW1);											// { i X_i V_i V_i^ V_i-V_i^ ���(S) STEP_i ���-�� ������� ���-�� �������� U_i |U_i - V_i| }

	double OLD_X = X;
	double OLD_V = V;
	double OLD_U = U;

	for (int i = 1; i <= parametrs.MAX_STEPS; ++i)
	{   
		double CURRENT_DOUBLING{};  // ���-�� �������� �� ������� ����
		double CURRENT_REDUCTION{}; // ���-�� ������� �� ������� ����
		
		double V_EXTRA = V;         // V^ ��� ����������� ����
		double X_EXTRA = X;			// X^ ��� ����������� ����

		for (int j = 0; j < 2; ++j) {
			make_Step(X_EXTRA, V_EXTRA, parametrs.STEP / 2.0);
		}
		make_Step(X, V, parametrs.STEP);

		double S = (V_EXTRA - V) / (pow(2.0, 4) - 1.0);
		ERRORS_LIST.emplace_back(abs(S * pow(2, 4)));
		
		double H = parametrs.STEP;
		int step_control = control_Error(X, V, X_EXTRA, V_EXTRA, OLD_X, OLD_V, S, CURRENT_DOUBLING, CURRENT_REDUCTION);

		double step;
		if (step_control == 0) {
			step = H;
		}
		else {
			step = parametrs.STEP;
		}
		double LP = abs(S * pow(2.0, 4));
		U = find_True_Solution(X, parametrs.START_POINT_FOR_U);

		if (X < parametrs.B - parametrs.E_BORDER) {
			std::vector<double> TABLE_ROW = { static_cast<double>(i), X, V, V_EXTRA, V - V_EXTRA, LP, step, CURRENT_REDUCTION, CURRENT_DOUBLING, U, abs(U - V) };
			TABLE_INFORMATION.emplace_back(TABLE_ROW);

			++reference.ITERATIONS_COUNT;
			STEPS_and_Xs.emplace_back(std::make_pair(step, X));
		}
		
		bool EXIT_FROM_FOR = false;
		if (x_in_border(parametrs.B,X,parametrs.E_BORDER)) {                       //�������� �� ��������� � ����������� ������ ������� �� X. 
			std::vector<double> TABLE_ROW = { static_cast<double>(i), X, V, V_EXTRA, V - V_EXTRA, LP, step, CURRENT_REDUCTION, CURRENT_DOUBLING, U, abs(U - V) };
			TABLE_INFORMATION.emplace_back(TABLE_ROW);

			++reference.ITERATIONS_COUNT;
			DISTANCE_Ui_Vi.emplace_back(abs(U - V));
			STEPS_and_Xs.emplace_back(std::make_pair(parametrs.STEP, X));

			EXIT_FROM_FOR = true;												   //���� X ����� � �����������, ��������� ��������������, ������ �� for �� �����
		}

		if (X > parametrs.B) {						// ���� X ����� �� ������ �������, ������������ �� ��� ����� � ������ ���,                                        
			X = OLD_X;								//  ������ ������� ������ ������� � ��������� ������ X � ����������� ��������������
			V = OLD_V;
			U = OLD_U;

			parametrs.STEP = parametrs.B - OLD_X;
			make_Step(X, V, parametrs.STEP);

			U = find_True_Solution(X, parametrs.START_POINT_FOR_U);
			DISTANCE_Ui_Vi.emplace_back(abs(U - V));

			std::vector<double> TABLE_ROW = { static_cast<double>(i), X, V, V_EXTRA, V - V_EXTRA, LP, parametrs.STEP, CURRENT_REDUCTION, CURRENT_DOUBLING, U, abs(U - V) };
			TABLE_INFORMATION.emplace_back(TABLE_ROW);
			++reference.ITERATIONS_COUNT;

			EXIT_FROM_FOR = true;                        //���� X ������ � ������ ��������, ����������� �������������� 
		}

		if (EXIT_FROM_FOR) {
			break;
		}
		else {
			DISTANCE_Ui_Vi.emplace_back(abs(U - V));
			OLD_X = X;
			OLD_V = V;
			OLD_U = U;
		}
	}

	reference.DISTANCE_B_LAST_POINT = parametrs.B - X;
	find_Max_Distance_Ui_Vi();
	find_Max_Step();
	find_Min_Step();
	find_Max_Error();
}

TestTask::FinalReference TestTask::get_reference()
{
	return reference;
}

std::list<std::vector<double>> TestTask::get_table_information()
{
	return TABLE_INFORMATION;
}

void TestTask::find_Max_Distance_Ui_Vi()
{ 
	reference.MAX_DISTANCE_U_V =  *std::max_element(DISTANCE_Ui_Vi.begin(), DISTANCE_Ui_Vi.end());
}





//������� ��� ��������
void TestTask::PrintTable(){
	if (TABLE_INFORMATION.front().size() == 6) {
		std::cout << std::setw(30) << "i" << std::setw(30) << "X" << std::setw(30) << "V" << std::setw(30) << "Step" << std::setw(30) << "U" << std::setw(30) << "|U-V|" << std::endl;
		for (auto it_list = TABLE_INFORMATION.begin(); it_list != TABLE_INFORMATION.end(); ++it_list) {
			for (size_t i = 0; i < 6; ++i) {
				if (i == 1) {
					std::cout << std::setw(30) << (*it_list)[i];
				}
				else {
					std::cout << std::setprecision(16) << std::setw(30) << (*it_list)[i];
				}
			}
			std::cout << std::endl;
		}
	}
	else {
		std::cout << std::setw(30) << "i" << std::setw(30) << "X" << std::setw(30) << "V" << std::setw(30) << "V^" << std::setw(30) << "V-V^" << std::setw(30) << "S" <<
			std::setw(30) << "Step" << std::setw(30) << "����������" << std::setw(30) << "��������" << std::setw(30) << "U_i" <<std::setw(30) << "|U_i-V_i|" << std::endl;
		for (auto it_list = TABLE_INFORMATION.begin(); it_list != TABLE_INFORMATION.end(); ++it_list) {
			for (size_t i = 0; i < 11; ++i) {
				if (i == 1) {
					std::cout << std::setw(30) << (*it_list)[i];
				}
				else if (i == 7 || i == 8) {
					std::cout << std::setw(30) << static_cast<int>((*it_list)[i]);
				}
				else {
					std::cout << std::setprecision(16) << std::setw(30) << (*it_list)[i];
				}
			}
			std::cout << std::endl;
		}
	}
}

void TestTask::Write_To_File()
{
	std::ofstream outFile("points.txt");
	if (!outFile) {
		std::cerr << "������ �������� �����!" << std::endl;
	}
	else {
		if (TABLE_INFORMATION.front().size() == 6) { //���� ������� ��� �������� ��
			for (auto it_list = TABLE_INFORMATION.begin(); it_list != TABLE_INFORMATION.end(); ++it_list) {
				outFile << std::fixed << std::setprecision(20) << (*it_list)[1] << "\t" << std::setw(30) 
						<< std::setprecision(20) << (*it_list)[2] << "\t" << std::setw(30) << std::setprecision(20) 
						<< (*it_list)[4] << std::endl;  // ������ � ������� "x v u"
			}
			outFile.close();
		}
		else { //���� ������� � ��������� ��
			for (auto it_list = TABLE_INFORMATION.begin(); it_list != TABLE_INFORMATION.end(); ++it_list) {
				outFile << std::fixed << std::setprecision(20) << (*it_list)[1] << "\t" << std::setw(30) 
						<< std::setprecision(20) << (*it_list)[2] << "\t" << std::setw(30) << std::setprecision(20) 
						<< (*it_list)[9] << std::endl;  // ������ � ������� "x v u"
			}
			outFile.close();
		}
	}
}

void TestTask::PrintReference()
{
	std::cout<< std::endl;
	std::cout << "DISTANCE_B_LAST_POINT :" << reference.DISTANCE_B_LAST_POINT << std::endl;
	std::cout << "ITERATIONS_COUNT :" << reference.ITERATIONS_COUNT << std::endl;
	std::cout << "MAX_DISTANCE_U_V :" << reference.MAX_DISTANCE_U_V << std::endl;
	std::cout << "MAX_ERROR :" << reference.MAX_ERROR << std::endl;
	std::cout << "MAX_STEP_WITH_X :" << "STEP = " << reference.MAX_STEP_WITH_X.first <<"  X = " << reference.MAX_STEP_WITH_X.second << std::endl;
	std::cout << "MIN_STEP_WITH_X :" << "STEP = " << reference.MIN_STEP_WITH_X.first <<"  X = " << reference.MIN_STEP_WITH_X.second << std::endl;
	std::cout << "STEP_DOUBLING_COUNT :" << reference.STEP_DOUBLING_COUNT << std::endl;
	std::cout << "STEP_REDUCTION_COUNT :" << reference.STEP_REDUCTION_COUNT << std::endl;
	std::cout << "IS_INF :" << reference.IS_INF << std::endl;
}





