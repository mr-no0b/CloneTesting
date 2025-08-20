#include<bits/stdc++.h>
using namespace std;

int main(){
    int x,y;cin>>x>>y;
    char c;cin>>c;
    if(c=='+'){
        cout<<x+y;
    }
    else if(c=='-'){
        cout<<x-y;
    }
    else if(c=='*'){
        cout<<x*y;
    }
    else if(c=='/'){
        cout<<x/y;
    }
}
