/*
 * Copyright (c) Valentin Galea
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
*/

/*
 * Largest integer constant available to the compiler
 *
 */
typedef unsigned long longest_t;

/*
 * The Class
 *
 */
template<longest_t N>
class bin
{
public:
        enum {
                value = (N % 8) + (bin<N / 8>::value << 1)
        };
};

/*
 * Specialization for ending the chain
 *
 */
template<>
class bin<0>
{
public:
        enum {
                value = 0
        };
};

/*
 * Macro-processing glue: force the number to be octal to both
 * end the recursion chain and make posible more digits
 *
 */
#define binary( n ) bin<0##n>::value

/*
 * Tests
 *
 */
/*void testing()
{
	int _3 = binary( 11 );
	int _7 = binary( 111 );
	int _15 = binary( 1111 );
	int _255 = binary( 11111111 );
	int _256 = binary( 100000000 );
	int _long_1 = binary( 0000000000001 );
}
*/
