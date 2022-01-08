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
struct Instance{    //实例
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
};
struct Tableins{    //表实例
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



void init_gridGraph(){
	gridGraph.resize(gridNum);
	for (int i = 0; i < gridGraph.size(); i++) {
		gridGraph[i].resize(gridNum);
	}
}


void getSpace(){
	for (int i = 0; i < instance.size(); i++) {
		int tempRow = instance[i].row / distance_pre;
		int tempCol = instance[i].col / distance_pre;
		gridGraph[tempRow][tempCol].push_back(instance[i]);
		instance[i].rowNum = tempRow;
		instance[i].colNum = tempCol;
	}
}
bool isNeighbor(double row1, double col1, double row2, double col2){
	double dis = sqrt(pow(row1-row2, 2) + pow(col1-col2, 2));
	if (dis <= distance_pre) {
		return true;
	}else {
		return false;
	}
}
void findNeibor_in_space(Instance& oriIns,  int rowNum, int colNum){
	for (int j = 0; j < gridGraph[rowNum][colNum].size(); j++) {
		if (isNeighbor(oriIns.row, oriIns.col, gridGraph[rowNum][colNum][j].row, gridGraph[rowNum][colNum][j].col)){
			if(oriIns.feature == gridGraph[rowNum][colNum][j].feature && oriIns.insNum == gridGraph[rowNum][colNum][j].insNum) continue;
			oriIns.true_neighbor.push_back(gridGraph[rowNum][colNum][j]);
			if (oriIns.feature < gridGraph[rowNum][colNum][j].feature) {
				oriIns.neighbor.push_back(gridGraph[rowNum][colNum][j]);
			}
		}
	}
}

void findNeighbor(){
	for (int i = 0; i < instance.size(); i++) {
		string strTemp = instance[i].feature + to_string(instance[i].insNum);
		instance[i].completeName = strTemp;
		int rowNum = instance[i].rowNum;
		int colNum = instance[i].colNum;
		findNeibor_in_space(instance[i], rowNum, colNum);
		if (rowNum != 0) {
			findNeibor_in_space(instance[i], rowNum-1, colNum);
		}
		if (rowNum != gridNum-1) {
			findNeibor_in_space(instance[i], rowNum+1, colNum);
		}
		if (colNum != 0) {
			findNeibor_in_space(instance[i], rowNum, colNum-1);
		}
		if (colNum != gridNum-1) {
			findNeibor_in_space(instance[i], rowNum, colNum+1);
		}
		if (rowNum != 0 && colNum != 0) {
			findNeibor_in_space(instance[i], rowNum-1, colNum-1);
		}
		if (rowNum != 0 && colNum != gridNum-1) {
			findNeibor_in_space(instance[i], rowNum-1, colNum + 1);
		}
		if (rowNum != gridNum-1 && colNum != 0) {
			findNeibor_in_space(instance[i], rowNum+1,colNum-1);
		}
		if (rowNum != gridNum-1 && colNum != gridNum-1) {
			findNeibor_in_space(instance[i], rowNum+1, colNum+1);
		}
	}
}

void gen_Table2ins(){
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
			}else {
				str_table[str].table.push_back(temp);
			}
			temp.clear();
		}
	}
}

void candi_PCP(){
	for (auto it = str_table.begin(); it != str_table.end(); it ++) {
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
		for(int i = 0;i < col_count.size(); i++){
			double PI =1.0 * col_count[i].size() / instance_sum[name[i]];
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
vector<string> partition2(vector<string>& head, vector<string>null_spcp2){
	vector<string> ans = head;
	set<string> temp;
	vector<string> true_ans;
	for (auto it = head.begin(); it != head.end(); it++) {
		for (int j = 0; j < null_spcp2.size(); j++) {
			char alph1 = null_spcp2[j][0];
			char alph2 = null_spcp2[j][1];
			if (it->find(alph1) != it->npos && it->find(alph2) != it->npos) {
				temp.insert(*it);
				string temp1 = *it, temp2 = *it;
				temp1.erase(it->find(alph1), 1);
				temp2.erase(it->find(alph2), 1);
				ans.push_back(temp1);
				ans.push_back(temp2);
			}
		}
	}
	for (int i = 0; i < ans.size(); i++) {
		bool flag = false;
		for (auto it = temp.begin(); it != temp.end(); it++) {
			if (ans[i] == *it) {
				flag = true;
			}
		}
		if (flag == false) {
			true_ans.push_back(ans[i]);
		}
	}
	return true_ans;
}
string gen_pp(string str){
	int max_spi = -1;
	string pp;
	for (int i = 0; i < str.size(); i++) {
		for (int j = i+1; j < str.size(); j++) {
			string cand_pp = "";
			cand_pp += str[i];
			cand_pp += str[j];
			if (str_table.count(cand_pp)) {
				if(str_table[cand_pp].PI > max_spi) {
					max_spi = str_table[cand_pp].PI;
					pp = cand_pp;
				}
			}
		}
	}
	return pp;
}
void gen_spcp3(string headL, string pp, set<string>& SPCP3){
	if (headL.size() == 3) {
		SPCP3.insert(headL);
		return;
	}
	string tempStr1 = headL;
	string tempStr2 = headL;
	tempStr1.erase(tempStr1.find(pp[0]),1);
	tempStr2.erase(tempStr2.find(pp[1]),1);
	string pp1 = gen_pp(tempStr1);
	string pp2 = gen_pp(tempStr2);
	gen_spcp3(tempStr1, pp1, SPCP3);
	gen_spcp3(tempStr2, pp2, SPCP3);
}

bool SNSI_isCover_C(Instance ins,set<char> spcp_feature){
	set<char> s;
	s.insert(ins.feature);
	for (int i = 0; i < ins.true_neighbor.size(); i++) {
		s.insert(ins.true_neighbor[i].feature);
	}
	for (auto it1 = spcp_feature.begin(); it1 != spcp_feature.end(); it1 ++) {
		bool flag = false;
		for (auto it2 = s.begin(); it2 != s.end(); it2 ++) {
			if (*it2 == *it1) {
				flag = true;
			}
		}
		
		if(flag == false) {
			return false;
		}
	}
	return true;
}

bool check_SPCP(string spcp){
	set<char> spcp_feature;
	for(auto c: spcp){
		spcp_feature.insert(c);
	}
	map<char, vector<string>> SPins;
	map<char, double> spr;
	for(auto c: spcp){
		for (int i = 0; i < instance.size(); i++) {
			if (instance[i].feature == c && SNSI_isCover_C(instance[i], spcp_feature)) {
				SPins[c].push_back(instance[i].completeName);
			}
		}
	}
	
	double min_spr = 10000;
	for(auto c: spcp){
		spr[c] = (1.0 * SPins[c].size())/instance_sum[c];
		if (spr[c] < min_spr) {
			min_spr = spr[c];
		}
	}
	if(min_spr >= min_pre){
		return true;
	}else {
		return false;
	}
}

bool is_can_merge(string& spcp1, string& spcp2){
	sort(spcp1.begin(), spcp1.end());
	sort(spcp2.begin(), spcp2.end());
	int count = 0;
	int i = 0;
	while (i != spcp1.size()) {
		if (spcp1[i] != spcp2[i]) {
			count ++;
		}
		i++;
	}
	if (count == 1) {
		return true;
	}else {
		return false;
	}
}
string mergeTwoStr(string str1, string str2){
	string ansStr;
	set<char> temp;
	for(auto c: str1){
		temp.insert(c);
	}
	for(auto c: str2){
		temp.insert(c);
	}
	for (set<char>::iterator it = temp.begin(); it != temp.end(); it ++) {
		ansStr += *it;
	}
	return ansStr;
}

void Merge(vector<string>& spcp, map <int, set<string> >& tree){
	bool flag = false;
	vector<string> next_spcp;
	if (spcp.size() <= 1) {
		return;
	}
	for(int i = 0; i < spcp.size(); i++){
		for(int j = i+1; j < spcp.size(); j++){
			if (is_can_merge(spcp[i], spcp[j])) {
				flag = true;
				string newSpcp =  mergeTwoStr(spcp[i], spcp[j]);
//				if(!check_SPCP(newSpcp)) continue;
				next_spcp.push_back(newSpcp);
				tree[newSpcp.size()].insert(newSpcp);
			}
		}
	}
	if(flag == false) return;
	if(next_spcp.size() > 1) Merge(next_spcp, tree);
}
void gen_MSPCP(map <int, set<string> > tree, int k, vector<string>& ans){
	set<string> candi_spcp = tree[k];
	for (auto it = candi_spcp.begin(); it != candi_spcp.end(); it ++) {
		if(check_SPCP(*it)){
			ans.push_back(*it);
		}
	}
}

void partition_based(){
	vector<string> head;
	vector<string> null_spcp2;
	set<string> SPCP3;
	map <int, set<string> > tree;
	string temp;
	for (auto it = pcp2.begin(); it != pcp2.end(); it ++) {
		string str = *it;
		char lastChar;
		if (it == pcp2.begin() || str[0] != lastChar) {
			head.push_back(temp);
			temp = str[0];
		}
		temp += str[1];
		lastChar = str[0];
	}
	head.erase(head.begin());
	//gen null_spcp2
	for (int i = 0; i < head.size(); i++) {
		string str = head[i];
		for (int j = 1; j < str.size(); j++) {
			int ascStr1 = str[j-1];
			int ascStr2 = str[j];
			if (ascStr1 != ascStr2 - 1) {
				char null_char = char(ascStr2-1);
				string null_str;
				null_str += str[0];
				null_str += null_char;
				null_spcp2.push_back(null_str);
			}
		}
	}
	vector<string> l = partition2(head,null_spcp2);
	for (int i = 0; i < l.size(); i++) {
		string pp = gen_pp(l[i]);
		gen_spcp3(l[i], pp, SPCP3);
	}
	vector<string> true_spcp3;
	for(auto it = SPCP3.begin(); it != SPCP3.end(); it ++){
		if (check_SPCP(*it)) {
			true_spcp3.push_back(*it);
			tree[it->size()].insert(*it);
		}
	}
	Merge(true_spcp3, tree);
	int maxSize;
	for (auto it = tree.begin(); it != tree.end(); it ++) {
		maxSize = it->first;
	}
	vector<string> ans;
	gen_MSPCP(tree, maxSize, ans);
	while (ans.size()==0) {
		gen_MSPCP(tree, maxSize--, ans);
	}
	for (int i = 0; i < ans.size(); i++) {
		cout << ans[i] << " ";
	}
//	for (auto it = tree.begin(); it != tree.end(); it ++) {
//		cout << "阶数:" << it->first << endl;
//		set<string> SPCP_set = it->second;
//		for(auto c: SPCP_set){
//			cout << c << " ";
//		}
//		cout << endl;
//	}
}
void gen_pcp2(string dataFile){
	ifstream myfile(dataFile);
	myfile >> min_pre;
	myfile >> distance_pre >> maxLength;
	Instance temp;
	set<string> ins_count;
	for (int i = 0; i < 10000; i++) {
		temp.insNum = -1;
		myfile >> temp.insNum >>temp.feature >> temp.row >> temp.col;
		if(temp.row > maxLength || temp.col > maxLength) continue;
		if (temp.insNum == -1) {
			break;
		}else {
			string strTemp = temp.feature + to_string(temp.insNum);
			ins_count.insert(strTemp);
			instance.push_back(temp);
		}
	}
	for (auto it = ins_count.begin(); it != ins_count.end(); it ++) {
		string str = *it;
		char feature = str[0];
		instance_sum[feature] ++;
	}
	gridNum = maxLength/distance_pre;
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
//	cout << insName_ins["B80"].neighbor.size();
//	for (int i = 0; i < insName_ins["B80"].neighbor.size(); i++) {   
//		cout << insName_ins["B80"].neighbor[i].completeName << " ";                        //存在无法输出邻居全名的问题
//	}
}
//void refresh(map<int, set<Instance> > & hashTable){
//	for (int i = 0; i < instance.size(); i++) {
//		hashTable[instance[i].degree].insert(instance[i]);
//	}
//}
bool cmp(Instance a, Instance b){
	if(a.degree < b.degree){
		return true;
	}else {
		return false;
	}
}
vector<Instance> getDeg(vector<Instance>& instance){
//	map<int, set<Instance> > hashTable;
//	vector<Instance> ans; 
//	for (int i = 0; i < instance.size(); i++) {
//		hashTable[instance[i].degree].insert(instance[i]);
//	}
//	int count = instance.size();
//	for (auto it = hashTable.begin(); it != hashTable.end(); it++) {
//		auto tempSet = it->second;
//		while(!tempSet.empty()){
//			auto inst = *(tempSet.begin());
//			tempSet.erase(tempSet.begin());
//			ans.push_back(inst);
//			for (int j = 0; j < inst.true_neighbor.size(); j++) {
//				Instance& to_dele_ins = inst.true_neighbor[j];
//				to_dele_ins.degree--;
//				refresh(hashTable);
//			}
//		}
//	}
//	return ans;
	sort(instance.begin(), instance.end(), cmp);
	return instance;
}

void MCHT(){
	vector<int>num_ins;
	const int len = instance.size();
	for (int i = 0; i < instance.size(); i++) {
		num_ins.push_back(i);
		instance[i].num_MCHT = i;
	}
	
	vector<Instance> deg_oder = getDeg(instance);
	vector<int> deg_num;
	for (int i = 0; i < deg_oder.size(); i++) {
		deg_num.push_back(deg_oder[i].num_MCHT);
	}
	
	for (int i = 0; i < deg_oder.size(); i++) {
		vector<bool> pre(len), last(len);
		int index = i;
		vector<bool> nebol_bool (deg_num.size());
		for (int j = 0; j < i; j++) {
			pre[deg_num[j]] = true;
			
		}
		for(int j = i+1; j < len; j++){
			last[deg_num[j]] = true;
		}
		for (int j = 0; j < instance[i].true_neighbor.size(); j++) {
			nebol_bool[instance[i].true_neighbor[j].num_MCHT] = true;
		}
		vector<bool>p(len), x(len);
		for(int j = 0; j < len; j++){
			p[j] = last[j]&nebol_bool[j];
			x[j] = pre[j]&nebol_bool[j];
		}
		
		
	}
	
}
int main(int argc, char *argv[]) {
	int begin = clock();
	gen_pcp2("/Users/mac/Desktop/data1.txt");
//	partition_based();
	
	int end = clock();
	cout << "程序耗时："  << 1.0*(end - begin)/CLOCKS_PER_SEC << "秒";
	
}