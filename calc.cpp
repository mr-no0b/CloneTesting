#include<iostream>
using namespace std;

int main(){
    int x,y;cin>>x>>y;
    cout<<"1.Add 2.Substract 3. Multiply 4.Devide"<<endl;
    int xx;cin>>xx;
    if(xx==1){
        cout<<x+y;
    }
    else if(xx==2){
        cout<<x-y;
    }
    else if(xx==3){
        cout<<x*y;
    }
    else if(xx==4){
        cout<<x/y;
        cout<<"hehe";
    }
}
//changed
