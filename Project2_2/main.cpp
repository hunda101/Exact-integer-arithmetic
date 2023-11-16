#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <bitset>
#include <cassert>
#include <vector>
#include <tuple>
using namespace std;
#include "big_int_function_def.h"

void print_menu();



class MultiplicationLong{

public:
    MultiplicationLong(){};
    
    bigint direct_multiplying(bigint num1, bigint num2){
        return num1*num2;
    }
    bigint fraction(bigint num, bigint divider){
        return num%divider;
    }
    bool is_complex_multiplication(bigint a, bigint b) {
        
        return a ==  bigint("0") || b ==  bigint("0") ? false : (a > bigint("10000") || b > bigint("10000"));
    }
    bigint shift(bigint num){
        bigint i ("0");
        bigint ten ("10");
        while (num > ten){
            
            num= num / ten;
            i+=1;
        }
        
        return i;
    }
    bigint find_max_pow(bigint num1, bigint num2) {
        
        bigint pow_1 = shift(num1);
        bigint pow_2 = shift(num2);
        
        bigint max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;
        return max_pow;
    }
    virtual bigint calculate(bigint num1, bigint num2){return 0;};
    bigint* multiplyPolynomials(bigint (*arr_arr)[3]) {
        bigint* result = new bigint[6]();
        for (int i = 0; i <= 4; ++i) {
            result[i] = 0;
        }

        for (int i = 0; i <= 2; ++i) {
            for (int j = 0; j <= 2; ++j) {
                result[i + j] += is_complex_multiplication(arr_arr[1][j], arr_arr[0][i]) ? calculate(arr_arr[1][j], arr_arr[0][i]) : (arr_arr[1][j] *  arr_arr[0][i]);
                
            }
            
        }
        return result;
    }

};
class Caracuba : public MultiplicationLong {
public:
    Caracuba(){};
    bigint calculate(bigint num1, bigint num2){
        
        bigint sum ("0");
        bigint arr [2] {num1, num2};
        bigint arr1[3] {};
        bigint arr2[3] {};
        bigint arr_arr[2][3] {arr1[3], arr2[3]};
        
        bigint max_pow = find_max_pow(num1, num2);
        
        bigint ten ("10");
        bigint x = bigint::_big_pow(ten, max_pow);
        find_common(arr, arr_arr, max_pow);
        bigint* ur = multiplyPolynomials(arr_arr);
        for (int i=0; i <= 3; ++i) {
            bigint big_i = bigint(i);
            sum += ur[i] *  bigint::_big_pow(x, big_i);
        }
        return sum;
    };
private:
    void find_common(bigint *arr, bigint (*arr_arr)[3], bigint max_pow) {
        for (int i =0; i<2; ++i) {
            bigint ten ("10");
            bigint right_side = bigint::_big_pow(ten, max_pow);
            bigint num1_fraction = fraction(arr[i], right_side);
            bigint left_side = arr[i]/right_side;
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = left_side;
            arr_arr[i][2] = 0;
         }
    }
    bigint shift(bigint num){
        bigint i = 0;
        bigint hundred ("100");
        cout << hundred << endl;
        while (num > hundred ){
            num/=10;
            i += 1;
        }
        return i;
    }
    
    
};
class ToomCook : public MultiplicationLong{
public:
    ToomCook(){};
    bigint calculate(bigint num1, bigint num2){
        
        bigint arr[2] {num1, num2};
        bigint arr1[3] {};
        bigint arr2[3] {};
        bigint arr_arr[2][3] {arr1[3], arr2[3]};
        bigint log_b_1 = big_log10(num1);
        bigint log_b_2 = big_log10(num2);
        bigint k = log_b_1/bigint(3) > log_b_2/bigint(3) ? log_b_1/bigint(3) + bigint(1) : log_b_2/bigint(3) +bigint(1);
        for (int i =0; i<2; ++i) {
            bigint ten("10");
            bigint second_pow = bigint::big_pow(ten, k) ;
            bigint num1_fraction = fraction(arr[i], second_pow);
            arr[i]/=second_pow;
            bigint second_side = arr[i]%second_pow;
            arr[i]/=second_pow;
            bigint first_side = arr[i];
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = second_side;
            arr_arr[i][2] = first_side;
        }

        return evaluating_values(arr_arr, k);
    };
    bigint evaluating_values(bigint (*arr_arr)[3], bigint max_pow){
        //evaluation
        bigint _p0 = arr_arr[0][0] +  arr_arr[0][2];
        bigint p0 = arr_arr[0][0];
        bigint p1 = _p0 + arr_arr[0][1];
        bigint p_1 = _p0 - arr_arr[0][1] ;
        bigint p_2 = ((p_1 +arr_arr[0][2])*bigint(2) -  arr_arr[0][0]);
        bigint p_inf = arr_arr[0][2];
        bigint _q0 = arr_arr[1][0] +  arr_arr[1][2];
        bigint q0 = arr_arr[1][0];
        bigint q1 = _q0 + arr_arr[1][1];
        bigint q_1 = _q0 - arr_arr[1][1] ;
        bigint q_2 = (q_1 +arr_arr[1][2])*bigint(2) -  arr_arr[1][0];
        bigint q_inf = arr_arr[1][2];

        //Pointwise multiplication
        bigint r0 = p0 * q0;
        bigint r1 = p1 * q1;
        bigint r_1 = p_1 * q_1;
        bigint r_2 = p_2 * q_2;
        bigint r_inf = p_inf * q_inf;
        //Interpolation
        bigint _r0 = r0;
        bigint _r4 = r_inf;
        bigint _r3 = (r_2 - r1)/bigint(3);
        bigint _r1 = (r1 - r_1)/bigint(2);
        bigint _r2 = (r_1 - r0);
        _r3 = (_r2 - _r3)/bigint(2) + bigint(2)*r_inf;
        _r2 = _r2 + _r1 - _r4;
        _r1 = _r1 - _r3;
        bigint *r_pointer1[5] = {&_r0, &_r1, &_r2, &_r3, &_r4};
        bigint sum("0");
        bigint ten("10");
        bigint power_of_ten = bigint::_big_pow(ten, max_pow);
        for(int i = 0; i < 5; ++i){
            bigint big_i = bigint(i);
            bigint x_i = bigint::_big_pow(power_of_ten, big_i);
            sum += *r_pointer1[i] * x_i;
        }
        return sum;
    }
    bigint find_max_pow(bigint num1, bigint num2) {
        bigint pow_1 = shift(num1);
        bigint pow_2 = shift(num2);
        bigint max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;

        return max_pow;

    }
};


class Strassen : public MultiplicationLong {
public:
    Strassen(){};
        void find_common(bigint *arr, bigint (*arr_arr)[3], bigint max_pow) {
        for (int i =0; i<2; ++i) {
            bigint ten("10");
            bigint two("2");
            bigint max_pow_doubled = max_pow*two;
            bigint second_pow = bigint::_big_pow(ten, max_pow);
            bigint first_pow = bigint::_big_pow(ten, max_pow_doubled);
            bigint num1_fraction = fraction(arr[i], second_pow);
            arr[i]/=second_pow;
            bigint second_side = arr[i];
            bigint first_side = arr[i] / first_pow;
            arr_arr[i][0] = num1_fraction;
            arr_arr[i][1] = second_side;
            arr_arr[i][2] = first_side;
        }
    }
    
    bigint calculate(bigint num1, bigint num2){
        
        bigint sum("0");
        bigint max_pow = find_max_pow(num1, num2);
        bigint arr [2] {num1, num2};
        bigint arr1[3] = {};
        bigint arr2[3] = {};
        bigint arr_arr[2][3] {arr1[3], arr2[3]};
        max_pow = !is_even(max_pow) ? max_pow+=1 : max_pow;
        bigint ten("10");
        bigint x = bigint::_big_pow(ten, max_pow);
        find_common(arr, arr_arr, max_pow);
        bigint* ur = multiplyPolynomials(arr_arr);
        for (int i=0; i <= 3; ++i) {
            bigint big_i = bigint(i);
            sum += ur[i] * bigint::_big_pow(x, big_i);
        }
        return sum;
    };
    bigint find_max_pow(bigint num1, bigint num2) {
        bigint pow_1 = min_pow_10(num1);
        bigint pow_2 = min_pow_10(num2);
        bigint max_pow = pow_1 >= pow_2 ? pow_1 : pow_2;
        return max_pow;
    }
    
    bool is_even(bigint num){
        
        return num % bigint('2') == bigint('0');
    }
    bigint min_pow_10(bigint num){
        bigint pow_10 = 0;
        bigint result = 1;
        bigint ten("10");
        while(result <= num){
            pow_10+=bigint(1);
            result = bigint::_big_pow(ten, pow_10);
            num/=bigint(10);
        }
        return pow_10;
    }
};


class Schonhage : public MultiplicationLong {
public:
    Schonhage(){};
    bigint calculate(bigint num1, bigint num2){
        
        bigint condition = 0;
        bigint m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, u1, u2, u3, u4, u5, u6,v1,v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6, uv1, uv2, uv3, uv4, uv5, uv6;
        for(int k = 0, j= 0 ; k < j; ++k){
            bigint q0 = 1;
            bigint q1 = 3*q0 -bigint("1");
            q0 = q1;
            m1 = bigint::_big_pow(bigint("2"), bigint("6")*q0-bigint("1"))-bigint("1");
            m2 = bigint::_big_pow(bigint("2"), bigint("6")*q0+bigint("1"))-bigint("1");
            m3 = bigint::_big_pow(bigint("2"), bigint("6")*q0+bigint("2"))-bigint("1");
            m4 = bigint::_big_pow(bigint("2"), bigint("6")*q0+bigint("3"))-bigint("1");
            m5 = bigint::_big_pow(bigint("2"), bigint("6")*q0+bigint("4"))-bigint("1");
            m6 = bigint::_big_pow(bigint("2"), bigint("6")*q0+bigint("5"))-bigint("1");
            bigint arr_arr [2] = {num1, num2};
            bigint* arr_m [6] = {&m1, &m2, &m3, &m4, &m5, &m6};
            for(int i = 0; i < 2; ++i){
                for(int j = 0; j < 6; ++j){
                    if(arr_arr[i] < *arr_m[j]){
                        cout << 1;
                        break;
                    }else{
                        condition = 1;
                    }
                    
                }
                if(condition == bigint(1)){
                    j+= 1;
                    break;
                }
            }
        }
        u1 = num1%m1;
        u2 = num1%m2;
        u3 = num1%m3;
        u4 = num1%m4;
        u5 = num1%m5;
        u6 = num1%m6;
        
        
        
        v1 = num2%m1;
        v2 = num2%m2;
        v3 = num2%m3;
        v4 = num2%m4;
        v5 = num2%m5;
        v6 = num2%m6;
        
        uv1 = v1*u1;
        uv2 = v2*u2;
        uv3 = v3*u3;
        uv4 = v4*u4;
        uv5 = v5*u5;
        uv6 = v6*u6;
        
        w1 = uv1%m1;
        w2 = uv2%m2;
        w3 = uv3%m3;
        w4 = uv4%m4;
        w5 = uv5%m5;
        w6 = uv6%m6;
        bigint freak = w4+m4;
        while (true){
            freak += m4;
            if(freak % m1 == w1 && freak % m2 == w2 && freak % m3 == w3 && freak % m4 == w4 && freak % m5 == w5 && freak % m6 == w6){
                break;
            }
        }
        return freak;
    }
    
    
    
    };
class BigDivision{
public:
    BigDivision(){};
protected:
    bigint to_binary(bigint t) {
        string binary_form = "";
        for(int i = 0; t != bigint("0"); ++i){
            binary_form = (t % bigint("2")).str + binary_form;
            t/=2;
        }
        return bigint::_to_bigint(binary_form);
    }
    string to_binary(int t) {
        string binary_form = "";
        for(int i = 0; t != 0; ++i){
            binary_form = to_string(t % 2) + binary_form;
            t/=2;
        }
        return binary_form;
    }
    
    string to_binary_counter(bigint t, bigint* counter) {
        string binary_form = "";
        for(int i = 0; t != bigint("0"); ++i){
            binary_form = (t % bigint("2")).str + binary_form;
            t/=2;
            *counter+=1;
        }
        
        return binary_form;
    }
    string move_left(bigint num){
        string moved_string = "0.";
        moved_string+= num.str;
        if(moved_string[2] != '1') moved_string[2] = 1;
        return moved_string;
    }
    string shift(string num, long long dot_pos){
        string moved_string = num;
        size_t dotPos1 = num.find('.');
        if(dotPos1 == std::string::npos){
            num = binary_floor(num, 0);
            dotPos1 = num.length()  - 1;
        }
        size_t newDotPos = dotPos1 >= dot_pos  ? dotPos1 - dot_pos : 0;
        moved_string.erase(dotPos1);
        moved_string.insert(newDotPos, ".");
        if (moved_string[0] == '.') {
            moved_string.insert(0, "0");
        }
        return (moved_string);
    }
    string binary_floor(string binary1, long long k){
        string subString = "";
        size_t dotPos1 = binary1.find('.');
        if(dotPos1 == std::string::npos){
            binary1 += ".";
            dotPos1 = binary1.length()  - 1;
            
        }
        size_t  len_of_decimal = binary1.length() - dotPos1 -1;
        if( dotPos1 <= len_of_decimal){
            binary1.erase(dotPos1 + k + 1);
            subString += binary1;
        }else{
            subString+= binary1;
            for(size_t i = 0; i < k - len_of_decimal; ++i){
                subString+="0";
            }
        }
        
        
        return subString;
    }
    
    
    
    void simplify_fraction(bigint &numerator, bigint &denominator) {
        bigint commonDivisor = bigint::gcd(numerator.str, denominator.str);
        numerator /= commonDivisor;
        denominator /= commonDivisor;
    }
    void simplify_fraction_to_denominator(bigint &numerator, bigint &denominator, bigint targetNumerator) {
        bigint commonDivisor = bigint::gcd(numerator.str, denominator.str);
        
        while(targetNumerator != denominator){
            denominator/=2;
            numerator/=2;
        }
    }
    string fraction_to_binary(bigint numerator, bigint denominator) {
        string binar = "";
        bigint intPart = numerator / denominator;
        if(intPart == bigint(0)) binar += "0";
        binar += to_binary(intPart).str;
        binar += ".";
        bigint remainder = numerator % denominator;
        string binar_remainder = "";
        while(remainder > bigint(0)){
            remainder *= 2;
            binar_remainder += (remainder / denominator).str;
            remainder %= denominator;
        }
        
        binar += binar_remainder;
        return binar;
    }
    bigint find_common_denominator(bigint a, bigint b) {
        bigint max = (a > b) ? a : b;
        while (true) {
            if (max % a == bigint(0) && max % b == bigint(0)) {
                return max;
            }
            max++;
        }
    }
    };

    class inversed_Cook: protected BigDivision{
    private:
        bigint num;
        
    public:
        inversed_Cook(bigint num) : num(num) {};
        tuple <string, bigint> calculate(){
            bigint counter(0);
            long long k = 0;
            string binar_fraction  = "";
            string binar_form = to_binary_counter(num.str, &counter);
            binar_form = move_left(binar_form);
            int v1 = binar_form[2] - '0', v2 = binar_form[3] - '0', v3 = binar_form[4] - '0';
            bigint z_nom = bigint(static_cast<int>((floor(32/(4*v1+2*v2+v3))/4)));
            bigint z_denom = 1;
            while(to_bigint(static_cast<int>(pow(2, k))) < counter){
                bigint vk_nom  =num ;
                bigint vk_denom = bigint::_big_pow(bigint(2), bigint(counter));
//                simplify_fraction(vk_nom, vk_denom);
                bigint z2_nom = z_nom*z_nom;
                bigint z2_denom = z_denom*z_denom;
//                simplify_fraction(z2_nom, z2_denom);
                
                bigint multiplication_nom = z2_nom * vk_nom;
                
                bigint multiplication_denom = vk_denom * z2_denom;
//                simplify_fraction(multiplication_nom, multiplication_denom);
                
                
                bigint common_denom = find_common_denominator(multiplication_denom, z_denom);
                bigint multiplier1 = common_denom / multiplication_denom;
                bigint multiplier2 = common_denom / z_denom;
                multiplication_nom *= multiplier1;
                multiplication_denom = common_denom;
                
                
                bigint z_times2_nom = z_nom * bigint(2)*multiplier2;
                bigint z_times2_denom = common_denom;
//                simplify_fraction_to_denominator(z_times2_nom, z_times2_denom, multiplication_denom);
                z_nom = z_times2_nom - multiplication_nom;
                z_denom = multiplication_denom;
//                simplify_fraction(z_nom, z_denom);
                binar_fraction =  fraction_to_binary(z_nom, z_denom);
                size_t dotPos1 = binar_fraction.find('.');
                for(size_t i =  dotPos1+1 + pow(2, k+1)   + 1; i < static_cast<long>(binar_fraction.length()); ++i){
                    if(binar_fraction[i-1] == '1'){
                        binar_fraction = binary_floor(binar_fraction, pow(2, k+1)+ 1);
                        binar_fraction[1 + pow(2, k+1)+ 1 ] = '1';
                        break;
                    };
                }
                k+=1;
                
                }
            return make_tuple(binar_fraction, counter);
        }
        
    };
    
    class Cook: public BigDivision {
    private:
        tuple <bigint, bigint> input_value(){
            bigint numerator(0), denominator(0);
            while(numerator == bigint(0)){
                cout << "Input numerator: ";
                cin >> numerator;
            }
            while(denominator == bigint(0)){
                cout << "Input denominator: ";
                cin >> denominator;
            }
            return make_tuple(numerator, denominator);
            
        }
        tuple<bigint, bigint> decimal_string_to_fraction(const string& decimalStr) {
            size_t decimalPointPos = decimalStr.find('.');
            string integerPart = decimalStr.substr(0, decimalPointPos);
            string decimalPart = decimalStr.substr(decimalPointPos + 1);
            bigint denominator = bigint::_big_pow(bigint(2), bigint::_to_bigint(static_cast<long>(decimalPart.length()))) ;
            bigint pow_of_two = denominator;
            bigint numerator = to_bigint("0");
            cout << endl;
            for(size_t i = 0; i < integerPart.length(); ++i){
                numerator += (integerPart[i]-'0')*pow_of_two;
                pow_of_two/=2;
            }
            for(size_t i = 0; i < decimalPart.length(); ++i){
                
                numerator += (decimalPart[i]-'0')*pow_of_two;
                pow_of_two/=2;
            }
            simplify_fraction(numerator, denominator);

            return make_tuple(numerator, denominator);
        }
        double binary_to_decimal(string binaryFraction) {
            int intPart = 0;
            double fracPart = 0.0;
            size_t dotPosition = binaryFraction.find('.');

            if (dotPosition != string::npos) {
                intPart = stoi(binaryFraction.substr(0, dotPosition), nullptr, 2);
            } else {
                intPart = stoi(binaryFraction, nullptr, 2);
            }

            if (dotPosition != string::npos && dotPosition + 1 < binaryFraction.length()) {
                string fracBinary = binaryFraction.substr(dotPosition + 1);
                for (size_t i = 0; i < fracBinary.length(); ++i) {
                    fracPart += (fracBinary[i] - '0') / pow(2, i + 1);
                }
            }

            return intPart + fracPart;
        }
        
    public:
        Cook(){};
        
        
        double calculate(){
            tuple<bigint, bigint> values_tup = input_value();
            bigint num1 = get<0>(values_tup);
            bigint num2 = get<1>(values_tup);
            inversed_Cook inversed_num(num2);
            tuple <string, bigint> inversedValue_tuple = inversed_num.calculate();
            string inversedValue = get<0>(inversedValue_tuple);
            tuple<bigint, bigint> fraction = decimal_string_to_fraction(inversedValue);
            bigint numerator = get<0>(fraction);
            bigint denominator = get<1>(fraction);
            numerator *= num1;
            denominator *= big_pow(bigint(2), get<1>(inversedValue_tuple));
            return binary_to_decimal(fraction_to_binary(numerator, denominator));
            
        }
        
        
    };
    
    
    
    class PrimeTestBase{
    public:
        PrimeTestBase(){};
        bool is_compared(long long d, long b, long long a){
            return ((d-b)%a==0);
        }
        bool odd(long long num){
            return num % 2 != 0;
        }
        long long gcd(long long a, long long b){
            long long t;
            while (b != 0){
                t = b;
                b = a%b;
                a = t;
            }
            return a;
        }
        long long jacobi(long long a, long long b){
            long long g;
            assert(odd(b));
            if(a >=b) a%=b;
            if(a == 0) return 0;
            if(a == 1) return 1;
            if(a < 0){
                if((b-1)/2 % 2 == 0){
                    return jacobi(-a, b);
                    
                }
                else{
                    return -jacobi(a, b);
                }
            }
            if(a % 2 == 0){
                if(((b*b-1)/8) % 2 ==0){
                    return +jacobi(a/2, b);
                }
                else{
                    return -jacobi(a/2, b);
                }
            }
            g = gcd(a, b);
            assert(odd(a));
            if(g == a){
                return 0;
            }
            else if(g != 1){
                return jacobi(g, b)* jacobi(a/g, b);
            }
            else if(((a-1)*(b-1)/4)%2 == 0){
                return + jacobi(b, a);
            }
            else{
                return - jacobi(b, a);
            }
        }
        bool is_power_two(long long num){
            for(int i = 1; i <= num; i*=2){
                if( i == num) return true;
            }
            return false;
        }
        
    };
    class Solovey_Strassen : public PrimeTestBase{
    public:
        Solovey_Strassen(){};
        double input(){
            long long num = 0, k;
            while (num == 0) {
                cout << "input odd number: ";
                cin >> num;
                if(!odd(num)){
                    num = 0;
                }
                
            }
            cout << "input accuracy: ";
            cin >> k;
            return calculate_probability(num, k);
            
        }
        double calculate_probability(long long p, long long k){
            for(long long i = 0; i< k; ++i){
                long long a = rand() % (p-3) + 2;
                if(gcd(a, p) > 1){
                    return 0;
                }
                long long j = static_cast<long long>(pow(a, (p-1)/2)) % p;
                long long Jacob = jacobi(a, p);
                long long r = (j - Jacob);
                r = r < 0 ? r += p : r;
                if((r % p) != p){
                    return 1;
                }
            }
            return 1 - pow(2, -k);
        }
        
    };
    
    class Luk_Lehmer : public PrimeTestBase{
    public:
        Luk_Lehmer(){};
        double input(){
            long long p = 0;
            while (p == 0) {
                cout << "input odd number and is power of two: ";
                cin >> p;
                if(!odd(p)){
                    p = 0;
                }
                if (!(is_power_two(p+1))){
                    p = 0;
                }
                
            }
            
            return calculate_probability(p);
            
        }
        double calculate_probability(long long p){
            long long S = 4;
            long long k = 1;
            long long M = static_cast<long long >(pow(2, p)) - 1;
            
            while(k != (p - 1)){
                S = ((S*S)- 2) % M;
                k += 1;
            }
            if(S == 0) return 1;
            return 0;
        }
        
    };
    class Rabin_Miller : public PrimeTestBase{
    public:
        Rabin_Miller(){};
        double input(){
            long long p = 0, k = 0;
            while (p == 0) {
                cout << "input odd number: ";
                cin >> p;
                if(!odd(p)){
                    p = 0;
                }
                
                
            }
            cout << "input accuracy: ";
            cin >> k;
            return calculate_probability(p, k);
            
        }
        double calculate_probability(long long p, long long k){
            
            int s = 0;
            long long m = p -1, t = m;
            while (t % 2 == 0) {
                t /= 2;
                s++;
            }
            long long d = m/pow(2, s);
            for (int j = 0; j < k; ++j) {
                long long a = rand() % (p-3) + 2;
                long long x = static_cast<long>(pow(a, d))%p;
                long long y= 0;
                for (int i = 0; i < s; ++i){
                    y = (x*x)%p;
                    if(y == 1 && x != 1 && x != (m)){
                        cout << "в первом" << endl;
                        return 0;
                    }
                    x = y;
                }
                if( y != 1) return 0;
            }
            return 1;
            
        }
    };
    int main() {
        print_menu();

            int type;
            cin >> type;
        
            switch (type) {
                case 1:{
                    bigint num1 (0), num2(0);
                    
                    while(num1 == bigint(0)){
                        cout << "Input first number: ";
                        cin >> num1;
                    }
                    while(num2 == bigint(0)){
                        cout << "Input second number: ";
                        cin >> num2;
                    }
                    Caracuba caracuba_method;
                    cout << caracuba_method.calculate(num1, num2) << endl;
                    break;
                }
                case 2:{
                    bigint num1 (0), num2(0);
                    
                    while(num1 == bigint(0)){
                        cout << "Input first number: ";
                        cin >> num1;
                    }
                    while(num2 == bigint(0)){
                        cout << "Input second number: ";
                        cin >> num2;
                    }
                    ToomCook toom_cook_method;
                    cout << toom_cook_method.calculate(num1, num2) << endl;
                    cout << toom_cook_method.direct_multiplying(num1, num2) << endl;
                    break;
                }
                case 3:
                {
                    bigint num1 (0), num2(0);
                    
                    while(num1 == bigint(0)){
                        cout << "Input first number: ";
                        cin >> num1;
                    }
                    while(num2 == bigint(0)){
                        cout << "Input second number: ";
                        cin >> num2;
                    }
                    Schonhage shonhage_method;
//                    cout << shonhage_method.calculate(a, b) << endl;
                    cout << shonhage_method.direct_multiplying(num1, num2) << endl;
                    break;
                }
                case 4:
                {
                    bigint num1 (0), num2(0);
                    
                    while(num1 == bigint(0)){
                        cout << "Input first number: ";
                        cin >> num1;
                    }
                    while(num2 == bigint(0)){
                        cout << "Input second number: ";
                        cin >> num2;
                    }
                    Strassen strassen_method;
//                    cout << strassen_method.calculate(a, b) << endl;
                    cout << strassen_method.direct_multiplying(num1, num2) << endl;
                    break;
                }
                case 5:
                {
                    bigint num(0);
                    
                    while(num == bigint(0)){
                        cout << "Input number: ";
                        cin >> num;
                    }
                    inversed_Cook inversed_cook(num);
                    cout << get<0>(inversed_cook.calculate()) << endl;
                    break;
                }
                case 6:
                {
                    Cook cook_method;
                    cout << cook_method.calculate() << endl;
                    break;
                }
                case 7:
                {
                    Luk_Lehmer luk_lehmer;
                    cout << luk_lehmer.input() << endl;
                    
                    break;
                }
                case 8:{
                    Rabin_Miller rabin_miller;
                    cout << rabin_miller.input() << endl;
                    break;
        
                }
                case 9:{
                    Solovey_Strassen solovey_strassen_method;
                    cout << solovey_strassen_method.input() << endl;
                    break;
                }
                default:{
                    cout << "Wrong input!!!\nOnly 1-9 requires";
                }
            }
        
            return 1;
        
    }
    
    
    void print_menu() {
        cout << "Enter number of task you want to handle with: " << endl;
        cout << "\t1 - Caratsuba multiplication of long numbers" << endl;
        cout << "\t2 - Toom-Cook multiplication of long numbers" << endl;
        cout << "\t3 - Shonhage multiplication of long numbers " << endl;
        cout << "\t4 - Strassen multiplication of long numbers" << endl;
        cout << "\t5 - Calculation inversed value with high accuracy" << endl;
        cout << "\t6 - Divising int numbers using Cook alghoritm" << endl;
        cout << "\t7 - Lehmer prime test" << endl;
        cout << "\t8 - Rabin-Miller prime test" << endl;
        cout << "\t9 - Solovey-Strassen prime test" << endl;
        
    }

