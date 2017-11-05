#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char const *argv[])
{
	vector<int> v;
	v.push_back(0);
	v.push_back(1);

	for (vector<int>::iterator i = v.begin(); i != v.end(); ++i)
	{
		cout << *i << endl;
	}
	return 0;
}