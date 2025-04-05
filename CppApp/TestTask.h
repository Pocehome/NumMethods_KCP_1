#pragma once
#include<list>
#include<vector>


class TestTask											// ����� ��������� �������� �����
{
private:
	void find_Max_Distance_Ui_Vi();						// ����� ��� ���������� ������������ ������� ����� Ui � Vi ��� �������� ������
    std::vector<double> DISTANCE_Ui_Vi;					// ������, ������� ������ ������� ����� Ui � Vi �� ������ ����
protected:
	struct MethodParametrs				// ��������� ��� ���������� ������
	{
		double A{};						// �������� �� ������� [a,b], �������� ����� ������� �� x. �������� �������������
		double B{};						// �������� �� ������� [a,b], �������� ������ ������� �� x. �������� �������������
		double STEP{};					// ��� ������
		double E_ERROR{};
		double E_BORDER{};				// ������� ��������� ��� ������� ������ �� ������� ������� [a,b]
		int MAX_STEPS{};				// ������������ ����� ����� ������
		double START_POINT_FOR_U{};		// ��������� ������� ��� U(x)
	};

	MethodParametrs parametrs;							// ������ ��������� MethodParametrs, � ������� ��������� �������� ������, ��������� � ��������� ����
	std::list<std::pair<double, double>> STEPS_and_Xs;	// ����, ������� ������ ���� � �������������� �� X
	std::list<double> ERRORS_LIST;						// ����, ������� ������ ��� �� ������ ����
	 
	bool x_in_border(const double& B, const double& X, const double& BORDER); //�������� �� ��������� X � ����������� ������ �������

	virtual double f(const double& X, const double& V);													        // ���������� ����� �������, � ������ ������ �� �������� ������ 
	double find_K1(const double& X, const double& V);													// �������, ������� ������� K1 �� ������
	double find_K2(const double& X, const double& V, const double& STEP, const double& K1);				// ...K2
	double find_K3(const double& X, const double& V, const double& STEP, const double& K2);				// ...K3
	double find_K4(const double& X, const double& V, const double& STEP, const double& K3);				// ...K4
	void make_Step(double& X, double& V, const double& STEP);                                           // �������, ����������� ���� ��� �������, ������ �������� X � V
	int control_Error(double& X, double& V, double& X_EXTRA, double& V_EXTRA, double& OLD_X,
					   double& OLD_V, double& S, double& CURRENT_DOUBLING, double& CURRENT_REDUCTION);  // �������, � ������� ���������� �������� ��, ������ ������ ���� �� ������

	void find_Max_Step();		// �������, ������� ���� ������������ ���, ������ ��������� MAX_STEP_WITH_X ��������� reference ������������� ���� (STEP,X)
	void find_Min_Step();       // ... ����������� ���
	void find_Max_Error();		// �������, ������� ���� ������������ ��, ������ ��������� MAX_ERROR ������������� ������������ ERROR �� ��������� parametrs
	void find_Min_Error();
	void find_Last_X();
	void find_Last_V();

public:
	struct FinalReference							 // ��������� ��� �������� �������
	{
		int ITERATIONS_COUNT{};						 // ������� ����� �������� 
		double DISTANCE_B_LAST_POINT{};				 // ���������� ����� ������ �������� �������������� B � ��������� ������ ��������������
		double MAX_ERROR{};							 // ������������ �������� ������ ���
		double MIN_ERROR{};							 // ����������� �������� ������ ���
		int STEP_DOUBLING_COUNT{};					 // ������� ����� �������� ����
		int STEP_REDUCTION_COUNT{};					 // ������� ����� ���������� ����
		std::pair<double, double> MAX_STEP_WITH_X{}; // ������������ ��� � �, ��� ������� ���� ��� ���������
		std::pair<double, double> MIN_STEP_WITH_X{}; // ����������� ��� � �, ��� ������� ���� ��� ���������
		double MAX_DISTANCE_U_V{};                   // ������������ ������� ����� U_i � V_i (������ ��� �������� ������)

		double LAST_X{};
		double LAST_V{};

		bool IS_INF{false};                          // ����, ������������ true , ���� ��� ���������� �������� ������������� �������� V
	};
	FinalReference reference;							// ������ ��������� Finalreference, � ������� ��������� �������� ������, ��������� � ��������� ����

	std::list<std::vector<double>> TABLE_INFORMATION;   // ������, ������� ������ ������ �������� �������
	
	TestTask(const double& _A, const double& _B, const double& _STEP, const double& _E_ERROR, 
			 const double& _E_BORDER, const int& _MAX_STEPS, const double& _START_POINT);    // ����������� ������
	~TestTask();																			 // ����������

	virtual void Solve_Without_Error_Control(); //����� ������� ��� ���
	virtual void Solve_With_Error_Control();	//����� ������� c ���

	FinalReference get_reference();                         //���������� ��������� reference
	std::list<std::vector<double>> get_table_information(); //���������� �������� �������

	// ������� ��� �������� � main , � �������� ��������, ��� ����� �� �����
	void virtual PrintTable();     //����� �������� ������� �� �������, ���� ����� ������ ������ �������, ����� �� �����������
	void virtual Write_To_File();  //������ X,Vi,Ui � ����
	void virtual PrintReference(); //����� �������� �������
};