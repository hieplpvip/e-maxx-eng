<!--?title Multiplicative functionn-->
# Multiplicative function

## Definition

In number theory, a multiplicative function is an arithmetic function $f(n)$ of a positive integer n that satisfies:

- $f(1) = 1$
- For every co-prime pair of integer $a$ and $b$: $f(ab)=f(a)f(b)$

Some popular multiplicative functions are:

- $Id_k(n)$: the power functions, defined by $id_k(n) = n^k$ for any number $k$. As special cases we have:
    + $Id_0(n) = I(n)$: the constant function
    + $Id_1(n) = Id(n)$: the [identity function](https://en.wikipedia.org/wiki/Identity_function)
- $gcd(n,k)$: the greatest common divisor of $n$ and a fixed integer $k$, as a function of $n$
- $\sigma_k(n)$: the divisor function, which is the sum of the k-th powers of all the positive divisors of $n$. As special cases we have:
    + $\sigma_0(n) = d(n)$: the number of positive divisors of $n$
    + $\sigma_1(n) = \sigma(n)$, the sum of all the positive divisors of $n$
- $\phi(n)$: [Euler's totient function](./algebra/phi-function.html)
- $\mu(n)$: the Mobius function
    + $\mu(n) = 0$ if $n$ is not square-free
    + $\mu(n) = 1$ if $n$ has even numbers of prime factors
    + $\mu(n) = -1$ if $n$ has odd numbers of prime factors

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

Here is the implementation for Euler's totient function:
```cpp
const int N = 10000000;
int lp[N + 1], phi[N + 1];
vector<int> pr;

phi[1] = 1;
for (int i = 2; i <= N; ++i) {
  if (lp[i] == 0) {
    // first case
    lp[i] = i;
    phi[i] = i - 1;
    pr.push_back(i);
  }
  for (int j = 0; j < (int)pr.size() && pr[j] <= lp[i] && i * pr[j] <= N; ++j) {
    lp[i * pr[j]] = pr[j];
    if (pr[j] < lp[i]) {
      // second case
      phi[i * pr[j]] = phi[i] * phi[pr[j]];
    } else {
      // third case
      phi[i * pr[j]] = phi[i] * pr[j];
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

## Applications
### Sum of GCD ### {#gcdsum}
Given a number $n$, find the sum of GCDs of all distinct pairs that can be formed with integers from $1$ to $n$.

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

## Practice Problems  
- [Codeforces - Bash Plays with Functions](https://codeforces.com/contest/757/problem/E)
- [Codeforces - The Holmes Children](https://codeforces.com/contest/776/problem/E)
