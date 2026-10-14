# Mathematical background

Let f: R^n -> R be differentiable. At a point w, the first-order change along
a direction d is grad f(w)^T d. Choosing d = -grad f(w) gives
-||grad f(w)||^2, which is strictly negative when the gradient is nonzero.
Thus a sufficiently small positive step decreases a smooth objective locally.

The iteration is

$$w_{k+1} = w_k - \alpha_k \nabla f(w_k).$$

For a function with an L-Lipschitz gradient, the descent lemma gives

$$f(w-\alpha\nabla f(w)) \le f(w)
 - \alpha(1-L\alpha/2)\|\nabla f(w)\|^2.$$

Hence 0 < alpha < 2/L ensures decrease at nonstationary points. This condition
requires information about the objective; MiniOpt does not estimate L or
promise convergence for arbitrary user-selected fixed steps.

For a symmetric positive definite matrix A,

$$f(w)=\tfrac12 w^T A w-b^T w,\qquad \nabla f(w)=Aw-b.$$

The unique minimizer solves Aw=b. With a fixed step between zero and
2/lambda_max(A), the gradient descent error contracts in the eigenbasis of A.
The examples compare the numerical result to `numpy.linalg.solve(A, b)`.

The diminishing schedule alpha_k = alpha_0/(k+1) reduces step lengths with
iteration count. Its behavior still depends on smoothness, scaling, and the
objective. It can converge slowly even on elementary quadratics.

Bounded numerical line search approximately minimizes f(w + alpha*d) on
[0, alpha_max] using SciPy's bounded scalar minimizer. Both endpoints are
compared with its candidate. On a nonconvex line, the candidate need not be a
global minimizer. A zero best step is a stopping event, not stationarity.

Autograd differentiates a composition of supported operations by the chain
rule. `autograd.numpy` supplies traceable operations. This is distinct from
central differences, which approximate a derivative with two nearby function
values and incur truncation and rounding errors.

MiniOpt declares success when ||grad f(w)|| <= tol. This is a stationarity
criterion: curvature or additional structure is required to classify the
point as a minimum. In particular, gradient descent can stop immediately
at a saddle point if initialized exactly there.
