#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <math.h>
#include <ctime>
#include <stack>
#include <unordered_map>
#include <algorithm>
using namespace std;
struct Instance {    //实例
	char feature;
	int insNum;
	string completeName;
	double row;
	double col;
	vector<Instance> neighbor;  //只计算字母序大于该特征的邻居
	vector<Instance> true_neighbor;
	int rowNum;     //格子行数
	int colNum;
	int degree;
	int num_MCHT;
	vector<Instance > childs;
};
struct Tableins {    //表实例
	string tableName;
	vector<vector<Instance> > table;
	int size;
	double PI;
};
double min_pre;   //参与度阈值
int distance_pre, maxLength;
int gridNum;

vector<Instance> instance;
vector<vector< vector<Instance>  > > gridGraph;
map<char, int> instance_sum;  //各特征的总数
map<string, Instance> insName_ins;
//set<char> all_feature;
map<string, Tableins> str_table;   //二阶表实例
set<string> pcp2;   //二阶频繁模式


void init_gridGraph() {
	gridGraph.resize(gridNum);
	for (int i = 0; i < gridGraph.size(); i++) {
		gridGraph[i].resize(gridNum);
	}
}


void getSpace() {
	for (int i = 0; i < instance.size(); i++) {
		int tempRow = instance[i].row / distance_pre;
		int tempCol = instance[i].col / distance_pre;
		gridGraph[tempRow][tempCol].push_back(instance[i]);
		instance[i].rowNum = tempRow;
		instance[i].colNum = tempCol;
	}
}
bool isNeighbor(double row1, double col1, double row2, double col2) {
	double dis = sqrt(pow(row1 - row2, 2) + pow(col1 - col2, 2));
	if (dis <= distance_pre) {
		return true;
	}
	else {
		return false;
	}
}
void findNeibor_in_space(Instance& oriIns, int rowNum, int colNum) {
	for (int j = 0; j < gridGraph[rowNum][colNum].size(); j++) {
		if (isNeighbor(oriIns.row, oriIns.col, gridGraph[rowNum][colNum][j].row, gridGraph[rowNum][colNum][j].col)) {
			if (oriIns.feature == gridGraph[rowNum][colNum][j].feature && oriIns.insNum == gridGraph[rowNum][colNum][j].insNum) continue;
			oriIns.true_neighbor.push_back(gridGraph[rowNum][colNum][j]);
			if (oriIns.feature < gridGraph[rowNum][colNum][j].feature) {
				oriIns.neighbor.push_back(gridGraph[rowNum][colNum][j]);
			}
		}
	}
}

void findNeighbor() {
	for (int i = 0; i < instance.size(); i++) {
		string strTemp = instance[i].feature + to_string(instance[i].insNum);
		instance[i].completeName = strTemp;
		int rowNum = instance[i].rowNum;
		int colNum = instance[i].colNum;
		findNeibor_in_space(instance[i], rowNum, colNum);
		if (rowNum != 0) {
			findNeibor_in_space(instance[i], rowNum - 1, colNum);
		}
		if (rowNum != gridNum - 1) {
			findNeibor_in_space(instance[i], rowNum + 1, colNum);
		}
		if (colNum != 0) {
			findNeibor_in_space(instance[i], rowNum, colNum - 1);
		}
		if (colNum != gridNum - 1) {
			findNeibor_in_space(instance[i], rowNum, colNum + 1);
		}
		if (rowNum != 0 && colNum != 0) {
			findNeibor_in_space(instance[i], rowNum - 1, colNum - 1);
		}
		if (rowNum != 0 && colNum != gridNum - 1) {
			findNeibor_in_space(instance[i], rowNum - 1, colNum + 1);
		}
		if (rowNum != gridNum - 1 && colNum != 0) {
			findNeibor_in_space(instance[i], rowNum + 1, colNum - 1);
		}
		if (rowNum != gridNum - 1 && colNum != gridNum - 1) {
			findNeibor_in_space(instance[i], rowNum + 1, colNum + 1);
		}
	}
}

void gen_Table2ins() {
	for (int i = 0; i < instance.size(); i++) {
		for (int j = 0; j < instance[i].neighbor.size(); j++) {
			string str;
			str += instance[i].feature;
			vector<Instance> temp;
			str += instance[i].neighbor[j].feature;
			temp.push_back(instance[i]);
			temp.push_back(instance[i].neighbor[j]);
			if (str_table.count(str) == 0) {
				Tableins tempTable;
				tempTable.table.push_back(temp);
				tempTable.size = str.size();
				tempTable.tableName = str;
				str_table[str] = tempTable;
			}
			else {
				str_table[str].table.push_back(temp);
			}
			temp.clear();
		}
	}
}

void candi_PCP() {
	for (auto it = str_table.begin(); it != str_table.end(); it++) {
		vector< vector<Instance> > table = it->second.table;
		string name = it->second.tableName;
		vector<set<string> > col_count;
		col_count.resize(it->second.size);
		for (int i = 0; i < table.size(); i++) {
			for (int j = 0; j < table[i].size(); j++) {
				string strTemp;
				strTemp += table[i][j].feature;
				strTemp += to_string(table[i][j].insNum);
				col_count[j].insert(strTemp);
			}
		}
		double minPI = 100000;
		for (int i = 0; i < col_count.size(); i++) {
			double PI = 1.0 * col_count[i].size() / instance_sum[name[i]];
			if (PI < minPI) {
				minPI = PI;
			}
		}
		str_table[name].PI = minPI;
		if (minPI >= min_pre) {
			pcp2.insert(name);
		}
	}
}
void gen_pcp2(string dataFile) {
	ifstream myfile(dataFile);
	myfile >> min_pre;
	myfile >> distance_pre >> maxLength;
	Instance temp;
	set<string> ins_count;
	for (int i = 0; i < 10000; i++) {
		temp.insNum = -1;
		myfile >> temp.insNum >> temp.feature >> temp.row >> temp.col;
		if (temp.row > maxLength || temp.col > maxLength) continue;
		if (temp.insNum == -1) {
			break;
		}
		else {
			string strTemp = temp.feature + to_string(temp.insNum);
			ins_count.insert(strTemp);
			instance.push_back(temp);
		}
	}
	for (auto it = ins_count.begin(); it != ins_count.end(); it++) {
		string str = *it;
		char feature = str[0];
		instance_sum[feature] ++;
	}
	gridNum = maxLength / distance_pre;
	init_gridGraph();
	getSpace();
	findNeighbor();
	gen_Table2ins();
	candi_PCP();
	for (int i = 0; i < instance.size(); i++) {
		string name = instance[i].completeName;
		insName_ins[name] = instance[i];
		instance[i].degree = instance[i].true_neighbor.size();
	}
}
Instance get_tree() {
	sort(instance.begin(), instance.end());
	Instance root;
	for (auto i : instance) {
		root.childs.push_back(i);
	}
	return root;
}
bool check(Instance& i, const set<Instance>& f) {
	for (auto ins : f) {
		if (find(i.true_neighbor.begin(), i.true_neighbor.end(), ins) == i.true_neighbor.end()) return false;
	}
	return true;
}
vector<set<Instance> > maxIns;
void gen_maximal(Instance& root, set<Instance>& f) {
	for (auto i : root.neighbor) {
		if (check(i, f)) {
			f.insert(i);
		}
	}
	maxIns.push_back(f);
}
map<string, int> count_candi;

int main()
{
	int co_minPre = 3;
	gen_pcp2();
	Instance root = get_tree();
	for (auto c : root.childs) {
		set<Instance> temp;
		gen_maximal(c, temp);
		temp.clear();
	}
	for (auto s : maxIns) {
		string tempStr = "";
		for (auto it = s.begin(); it != s.end(); it++) {
			tempStr += it->feature;
		}
		sort(tempStr.begin(), tempStr.end());
		count_candi[tempStr] ++;
	}
	for (auto it = count_candi.begin(); it != count_candi.end(); it++) {
		if (it->second >= co_minPre) {
			cout << it->first << endl;
		}
	}
}
