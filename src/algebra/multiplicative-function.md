<!--?title Multiplicative functionn-->
# Multiplicative function

## Definition

In number theory, a multiplicative function is an arithmetic function $f(n)$ of a positive integer n that satisfies:

- $f(1) = 1$
- For every co-prime pair of integer $a$ and $b$: $f(ab)=f(a)f(b)$

### Examples

- $Id_k(n)$: the power functions, defined by $id_k(n) = n^k$ for any number $k$. As special cases we have:
    + $Id_0(n) = I(n)$: the constant function
    + $Id_1(n) = Id(n)$: the [identity function](https://en.wikipedia.org/wiki/Identity_function)
- $\epsilon(n)$: the unit function, defined by $\epsilon(n) = 1$ if $n = 1$ and $0$ otherwise.
- $gcd(n,k)$: the greatest common divisor of $n$ and a fixed integer $k$, as a function of $n$
- $\sigma_k(n)$: the divisor function, which is the sum of the k-th powers of all the positive divisors of $n$. As special cases we have:
    + $\sigma_0(n) = d(n)$: the number of positive divisors of $n$
    + $\sigma_1(n) = \sigma(n)$, the sum of all the positive divisors of $n$
- $\phi(n)$: [Euler's totient function](./algebra/phi-function.html)
- $\mu(n)$: the Mobius function
    + $\mu(n) = 0$ if $n$ is not square-free
    + $\mu(n) = 1$ if $n$ has even numbers of prime factors
    + $\mu(n) = -1$ if $n$ has odd numbers of prime factors

### Transitive property:

Given two multiplicative functions $f(x)$ and $g(x)$, we have:

- Their product $h(x) = f(x) g(x)$ is also multiplicative.
- Their Dirichlet convolution $$(f * g)(n) = \sum_{d | n}{f(d) * g(n/d)}$$ is also multiplicative (more about Dirichlet convolution below).  
Specifically, if $f(x)$ is multiplicative, $h(x) = \displaystyle \sum_{d | n}{f(d)}$ is also multiplicative.

## Calculation using the sieve of Eratosthenes

Let's consider this problem:

- Given $n \le 1e7$ and a function $f(x)$, calculate the value of $f(i)$ for every integer $i$ in the segment $[1,n]$.

If we can prove that $f(x)$ is multiplicative and can calculate $f(p^k)$ in $O(1)$ for any prime number $p$ and integer $k$, we can utilize the [linear prime sieve](./algebra/prime-sieve-linear.html) to solve this problem in $O(n)$.

Consider the original implementation of linear prime sieve:
```cpp
const int N = 10000000;
int lp[N + 1];
vector<int> pr;

for (int i = 2; i <= N; ++i) {
  if (lp[i] == 0) {
    lp[i] = i;
    pr.push_back(i);
  }
  for (int j = 0; j < (int)pr.size() && pr[j] <= lp[i] && i * pr[j] <= N; ++j)
    lp[i * pr[j]] = pr[j];
}
```

For simplicity, let denote $pr[j]$ as $p$. There are $3$ cases:

- $i$ is prime ($lp[i] = 0$). In this case, we can determine the value of $f(i)$ directly.
- $p < lp[i]$. Since $lp[i]$ is the smallest prime factor of $i$, $i$ and $p$ must be co-prime.  
We can assign $f(ip) = f(i)f(p)$
- $p = lp[i]$. This is the most complicated case since $i$ and $p$ are not co-prime.  
In most situation, a simple relationship between them exists. For example, with Euler's totient function, we can easily infer $\phi(ip) = \phi(i) * p$ using the [fourth property](./algebra/phi-function.html#toc-tgt-0).

**Note**: Don't forget to initialize $f[1] = 1$ (very important).

Here is the to calculate both Euler's totient function and Mobius function:
```cpp
const int N = 10000000;
int lp[N + 1], phi[N + 1], mu[N + 1];
vector<int> pr;

phi[1] = 1;
for (int i = 2; i <= N; ++i) {
  if (lp[i] == 0) {
    // first case
    lp[i] = i;
    phi[i] = i - 1;
    mu[i] = -1;
    pr.push_back(i);
  }
  for (int j = 0; j < (int)pr.size() && pr[j] <= lp[i] && i * pr[j] <= N; ++j) {
    lp[i * pr[j]] = pr[j];
    if (pr[j] < lp[i]) {
      // second case
      phi[i * pr[j]] = phi[i] * phi[pr[j]];
      mu[i * pr[j]] = mu[i] * mu[pr[j]];
    } else {
      // third case
      phi[i * pr[j]] = phi[i] * pr[j];
      mu[i * pr[j]] = 0;
    }
  }
}
```

In case we can't find any relationship between $i$ and $p$, we can workaround by using an auxiliary array $cnt[i]$ which denotes the power of $p$ in $i$.

As $\dfrac{i}{p^{cnt[i]}}$ and $p^{cnt[i] + 1}$ are co-prime, we have:
$$f(ip) = f\left (\dfrac{i}{p^{cnt[i]}}  \right ) * f\left (p^{cnt[i] + 1}  \right )$$

To avoid unnecessarry recalculation of $p^{cnt[i]}$, we will store it as $pw[i]$.

Refer to [sum of GCD application](#gcdsum) below for sample implementation.

## Dirichlet convolution

[Dirichlet convolution](https://en.wikipedia.org/wiki/Dirichlet_convolution) is a way to generate a new function from two functions.

Formally, the Dirichlet convolution of two functions $f(x)$ and $g(x)$ is:
$$(f * g)(n) = \sum_{d_1 * d_2=n}{f(d_1) * g(d_2)}$$

Another representation is:
$$(f * g)(n) = \sum_{d | n}{f(d) * g(n/d)}$$

It has a very nice property: the Dirichlet convolution of two multiplicative functions is also a multiplicative function.
<details>
  <summary>Proof</summary>
  
  Consider any two co-prime numbers $a$ and $b$. Every divisor $d$ of $ab$ can be uniquely reprensented as $rs$, with $r|a$, $s|b$ and $gcd(r,s) = 1$.
  
  $$
  \begin{align}
  (f * g)(ab) &= \sum_{r|a, s|b}{f(rs) g(ab/rs)} \\\\
  &= \sum_{r|a, s|b}{f(r) f(s) g(a/r) g(b/s)} \\\\
  &= \sum_{r|a}{f(r) g(a/r)} \sum_{s|b}{f(s) g(b/s)} \\\\
  &= (f * g)(a) (f * g)(b)
  \end{align}
  $$
  
  Therefore, $(f * g)$ is multiplicative.
</details>

### Example
To illustrate the power of Dirichlet convolution, we will use it to prove the divisor function is multiplicative.  
Let $f(n) = n^k$ and $g(n) = 1$. Obviously $f$ and $g$ are multiplicative functions.
$$(f * g)(n) = \sum_{d | n}{f(d) * g(n/d)} = \sum_{d | n}{d^k} = \sigma_k(n)$$
Therefore, the divisor function is multiplicative.

You can find more examples on [Wikipedia](https://en.wikipedia.org/wiki/Multiplicative_function#Convolution).

## Applications
### Sum of GCD ### {#gcdsum}
Given a number $n$, find the sum of GCDs of all unordered pairs that can be formed with integers from $1$ to $n$ (unordered means that $(i, j)$ and $(j, i)$ are considered the same).

We will calculate $h(x) = \displaystyle \sum_{1 \le i \le x}{gcd(x, i)}$.  
The answer for the problem would be $\displaystyle \sum_{1 \le i \le n}{\left (h(i) - i  \right )}$

$h(x)$ can be rewritten as: $h(x) = \displaystyle \sum_{d | x}{d * count(d)} = \sum_{d | x}{d * \phi(x / d)}$

Here $count(d)$ is the number of pairs $(x, i)$ with GCD equal to $d$. For every such pair, we have: $gcd(x / d, i / d) = 1$.  
Therefore, $count(d) = \phi(x / d)$.

Let's use Dirichlet convolution with $f(x) = x$ and $g(x) = \phi(x)$. We get $h(x)$ is multiplicative.

The only thing left now is how to calculate $h(p^k)$ for any prime number $p$ and integer $k$.

$$h(p^k) = \sum_{d | p^k}{d * \phi(p^k / d)} = p^k + \sum_{0 \le i \le k-1}{p^i(p^{k-i}-p^{k-i-1})} = p^k + k * (p^k - p^{k - 1}) = (k + 1) * p^k - k * p^{k - 1}$$

Sample code:
```cpp
const int N = 10000000;

int lp[N + 1], cnt[N + 1], pw[N + 1];
long long h[N + 1], ans[N + 1];
vector<int> pr;

h[1] = 1; ans[1] = 0;
for (int i = 2; i <= N; ++i) {
  if (lp[i] == 0) {
    // first case
    lp[i] = pw[i] = i;
    cnt[i] = 1;
    h[i] = 2 * i - 1;
    pr.push_back(i);
  }
  ans[i] = ans[i - 1] + h[i] - i;
  for (int j = 0; j < (int)pr.size() && pr[j] <= lp[i] && i * pr[j] <= N; ++j) {
    lp[i * pr[j]] = pr[j];
    if (pr[j] < lp[i]) {
      // second case
      cnt[i * pr[j]] = 1;
      pw[i * pr[j]] = pr[j];
      h[i * pr[j]] = h[i] * h[pr[j]];
    } else {
      // third case
      cnt[i * pr[j]] = cnt[i] + 1;
      pw[i * pr[j]] = pw[i] * pr[j];
      // tmp = h[pr[j] ^ (cnt[i] + 1)]
      ll tmp = 1LL * (cnt[i * pr[j]] + 1) * pw[i * pr[j]] - 1LL * cnt[i * pr[j]] * pw[i];
      h[i * pr[j]] = h[i / pw[i]] * tmp;
    }
  }
}
```

### Number of co-prime pairs in a set

Given a set $S$ of numbers. Find the number of co-prime pairs in this set.

We will use inclusion-exclusion principle and Mobius function to solve this problem.

For simplicity, denote a pair $(x, y)$ divisible by $k$ if both $x$ and $y$ are divisible by $k$.

Using inclusion-exclusion principle, the answer is:
```
+ number of pairs divisible by 1
- number of pairs divisible by 2
- number of pairs divisible by 3
- number of pairs divisible by 5
...
+ number of pairs divisible by 2 * 3
+ number of pairs divisible by 2 * 5
+ number of pairs divisible by 3 * 5
...
- number of pairs divisible by 2 * 3 * 5
...
```

Let $cnt_k$ be the number of multiples of $k$ in $S$.

The answer is: $\displaystyle \sum{\dfrac{cnt_i * (cnt_i - 1)}{2} * \mu(i)}$.

Some problems require adding/removing numbers in a set. To add/remove $x$, we only need to update $cnt_k$ for $k$ divides $x$. Complexity of each operation is $O(\sigma_0(x))$.

#### Arbitrary GCD

To count pairs $(x, y)$ with $gcd(x, y) = g$ for any $g$, take all multiple of $g$ in $S$ and divide them by $g$ to create a new set. We move back to the original problem: count co-prime pairs in the new set.

#### Number of relative primes

Instead of counting pairs, we want to know how many numbers in $S$ are relative prime to $x$.

Still using inclusion-exclusion principle and Mobius function, the formula is:
$$\sum_{d | x}{cnt_d * \mu(d)}$$

#### Optimize with bitmask

When we update $cnt$ for some number $x$, we consider many unnecessary divisors of $x$: non square-free divisors, for which Mobius function is $0$.

For example, with $x = 16 = 2^4$, we have 5 divisors ($1$, $2$, $4$, $8$, $16$), but only two of them count ($1$ and $2$). The rest have Mobius value of $0$.

Let the prime factorization of $x$ be $p_1^{e_1} \cdot p_2^{e_2} \cdots p_k^{e_k}$. We will only consider divisors of type $p_1^{a_1} \cdot p_2^{a_2} \cdots p_k^{a_k}$ where $a_i \in \left \\{ 0; 1 \right \\}$. It's clear that we can use bitmask, each bit corresponds to one prime factor.

But how many bits will we need to consider? The number with most prime factors is $2 * 3 * 7 * \cdots$. Keep multiplying until it's larger than maximum of $x$.

For example, if maximum of $x$ is $5e5$, we only need to consider at most $6$ bits, since $2 * 3 * 5 * 7 * 11 * 13 * 17 = 510510$ is larger than $5e5$.

In addition, we don't need to pre-calculate Mobius function. Count number of bits set, if it's odd, Mobius value is $-1$, otherwise it's $1$.

## Practice Problems  
- [Codeforces - Bash Plays with Functions](https://codeforces.com/contest/757/problem/E)
- [Codeforces - The Holmes Children](https://codeforces.com/contest/776/problem/E)
- [Codeforces - Mike and Foam](https://codeforces.com/contest/547/problem/C)
- [Codeforces - Classical?](https://codeforces.com/contest/1285/problem/F)
- [HackerRank - Pairwise GCD](https://www.hackerrank.com/contests/infinitum13/challenges/pairwise-gcd)
