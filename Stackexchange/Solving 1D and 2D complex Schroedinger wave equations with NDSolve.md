**Solving 1D and 2D complex Schroedinger wave equations with NDSolve**

I do not agree with you when you write:

> I know the NDSolve is not magic...

My opinion is that *NDSolve* is one of the most complex functionality I've met so far in the Mathematica environment, with its millions of options and special function this is a real complex thing and it is hard indeed to get a proper result. NDSolve is big, you could even sell this as a standalone program.

This post is somewhat long and at one side it is an answer (or an approximation of answer) and on the other side it is opening a huge amount of new questions...

It's a pity that guys like *Rob Knapp* are not available to ask. So it's only me who tries to give a look under the hood of NDSolve while scratching only at the surface.)


----------


*NDSolve* includes a general solver for partial differential equations based on the methods of lines. There are several ways to control the selection of the spatial grid (I'll list here only those who are relevant for the 1D case):

> *AccuracyGoal          $\rightarrow$          the number of digits of absolute tolerance*

> *PrecisionGoal         $\rightarrow$          the number of digits of relative tolerance*

> *MinStepSize           $\rightarrow$          the minimum grid spacing to use*

> *MaxStepSize           $\rightarrow$          the maximum grid spacing to use*

The discretization is done (not for pseudospectral methods) with uniform grids with the following direct correspondences (with interval length L):

> MaxPoints $\rightarrow$ n    $\Leftrightarrow$    MaxStepSize $\rightarrow$ L/n

> MinPoints $\rightarrow$ n    $\Leftrightarrow$    MinStepSize $\rightarrow$ L/n

**1) The 1D case**:

    sol = NDSolve[{I D[u[t, x], t] == -D[u[t, x], {x, 2}],
        u[0., x] == Exp[-(x^2.)], u[t, 5.] == 0, u[t, -5.] == 0},
        u, {t, 0., 20.}, {x, -5., 5.}, MaxStepSize -> 0.01,
        AccuracyGoal -> 3, PrecisionGoal -> 3]

    Animate[Plot[Evaluate[Abs[u[t, x] /. First[sol]]^2], {x, -5, 5},
        PlotRange -> {0, 1}], {t, 0, 17, 0.01}]

![enter image description here][1]

In order to increase the grid size you have to decrease *MaxStepSize*.
Now your plot will work nicely.


----------


**Warming up**:

Now let's turn to a 2D heat equation. To get warm we'll solve a SIAM 100 Challenge using *NDSolve* for that, although an analytical solution would be the best way.

(I only want to show that the solver is quite capable to get an answer)

The SIAM 100 Challenge #8 states the following problem:

> A square plate [-1,1] is at temperature u = 0. At time t = 0 the
temperature is increased to u = 5 along one of the four sides while being
held at u = 0 along the other three sides, and heat then flows into the plate
according $u_i = \Delta{u}$. When does the temperture reach u = 1 at the center
of the plate? (Folkmar Bornemann)

The initial condition in that problem is discontinuous, so we'll provide a specific grid spacing.

    Quiet[Block[{n = 25, he = 5 UnitStep[-(x + 1)]}, hsol = NDSolve[
    {
     D[u[t, x, y], t] == D[u[t, x, y], x, x] + D[u[t, x, y], y, y],
     he == u[0, x, y],
     u[t, -1, y] == 5, u[t, 1, y] == 0, u[t, x, -1] == he,
     u[t, x, 1] == he
     }, u, {t, 0, 1}, {x, -1, 1}, {y, -1, 1},
     Method -> {"MethodOfLines",
      Method -> {"EventLocator",
        "Event" -> u[t, 0, 0] - 1,
        "EventAction" :> Throw[end = t, "StopIntegration"]},
        "SpatialDiscretization" -> {"TensorProductGrid",
            "MinPoints" -> {n, n}, "MaxPoints" -> {n, n}}}]]]

which yields an answer quickly:

    u->InterpolatingFunction[{{0.,0.424014},{-1.,1.},{-1.,1.}},<>]}}

We've specified the *EventLocator* controller method. Every time the *Event* option is zero *EventAction* is evaluated, which in this case will stop the integration.

The most important suboptions for method option = *MethodOfLines* are:

> "SpatialDiscretization $\rightarrow$ "TensorProductGrid"

> "MinPoints" $\rightarrow$ list Of Discreatization points for each spatial variable

> "MaxPoints" $\rightarrow$ dito

> "DifferenceOrder" $\rightarrow$ positive integer or "PseudoSpectral"


Now let's *Plot3D* the solution of the heat equation at t = 0.424014, along with a *DensityPlot*:

    GraphicsGrid[
    {{
        Plot3D[Evaluate[u[0.424014, x, y] /. hsol], {x, -1, 1}, {y, -1, 1}],
        DensityPlot[Evaluate[Abs[u[0.424014, x, y]] /. hsol], {x, -1, 1}, {y, -1, 1},
            PlotPoints -> 200, Mesh -> False]
    }}
    ]

![enter image description here][2]


----------


**The 1D complex Schroedinger wave equation**

We will model the quantum-mechanical scattering of a Gaussian wave packet by using a 1D time-dependent Schroedinger equation.

Let's define the the Schroedinger equation:

The time-dependent Schroedinger equation is defined as:

> $i\hbar \frac{d}{dt} \Psi(x,t)=[\frac{-\hbar^2}{2m} \Delta^2+V(x,t)] \psi(x,t)$

As you already wrote in your *Note 1*, the inital condition is a simple spreading Gaussian probability.

And its potential/kinetic energy is defined as:

> $e^{-(x+3)^{2}+3 i x}$

If we integrate this with $\Psi$ we get the maximum value of the potential, which is $\approx 6$

Setting $\hbar$ and $m$ = 1, and the potential to be 6, we get:

    schroedingerEq = I D[u[x, t], {t, 1}] ==  -1/2 D[u[x, t], {x, 2}] + 6
        Exp[-x^2] u[x, t]

In order to make sure that at $\pm xMax$ our wave function behaves identically we define a Dirichlet boundary condition which adds a cosine to `u`:

    DirichletBC[u_, x_, xM_] := u - (((u /. x -> -xM) - (u /. x -> xM))/
      2*(Cos[(x + xM)/(2 xM) Pi] + 1) + (u /. x -> xM))

Now let's solve this wave equation numerically:

    With[{xMax = 15},(nsol =
        NDSolve[{schroedingerEq,
            u[x, 0] == DirichletBC[Exp[-(x + 3)^2], x, xMax] Exp[3 I x],
            u[xMax, t] == 0, u[-xMax, t] == 0},
            u[x, t], {x, -xMax, xMax}, {t, 0, 5},
            AccuracyGoal -> 3, PrecisionGoal -> 3]) // Timing]

Plotting `nsol` yields:

    DensityPlot[Evaluate[Abs[u[x, t]] /. nsol], {x, -15, 15}, {t, 0, 5},
        PlotPoints -> 200, Mesh -> False]

Here we can see that at t $\approx$ 3.2 the wave packet reaches the right boundary and gets reflected there.

![enter image description here][3]


----------


**The 2D complex Schroedinger wave equation**

When we look at the byte size of `sol` (see the 1D part of this answer) we realize, that the size is about 72 MB large.

In order to avoid such huge *InterpolatingFunction* objects and since I found out, in our prior discussion, that a low memory consumption is important to you I want to show another way to tackle the problem.

Under the hood *NDSolve* is broken up actually in several parts:

 1. Equation processing and method selection
 2. Method initialization
 3. Numerical solution
 4. Solution processing

Normally, if you use *NDSolve* you won't even notice these steps, but there exist low-level functions which you can use on your own, to break up the steps.

 - ProcessEquations
 - Iterate
 - ProcessSolutions

*ProcessEquations* is used to set up the problem and to create the *StateData* data structure.

*Iterate* advances the numerical solution.

*ProcessSolutions* converts the numerical data into a *InterpolationFunction*

[Please see the documentation for NDSolve steps and components][4]

Let's first set up the problem using a homogenous Dirichlet boundary condition on all four edges, with length `L`, the initial condition `ics0` and the DifferenceOrder option `do`:

    makeNDSolveStateDataObject[ics0_, do_, L_: 5, opts___] :=
        First[NDSolve`ProcessEquations[{
        I D[u[t, x, y], t] == -D[u[t, x, y], {x, 2}] - D[u[t, x, y], {y, 2}],
        u[0., x, y] == ics0,
        u[t, -L, y] == u[t, L, y], u[t, x, -L] == u[t, x, L]},
        u, {t, 0., 2.}, {x, -L, L}, {y, -L, L}, opts,
        Method -> {"MethodOfLines",
             "SpatialDiscretization" -> {"TensorProductGrid",
             "DifferenceOrder" -> do}}]]

and build the *StateData* data structure:

    stateD = makeNDSolveStateDataObject[Exp[-(x^2 + y^2)], "Pseudospectral"];

    stateD ==> NDSolve`StateData[<0.>]

The nice thing about the *StateData* data structure is that we have more control of the integration. Sometimes it is appropriate to check the solution and change maybe some parameters and to start it over again.

*Iterate* on its own does not return a value but it modifies the *StateData* data structure.

We can integrate using *Iterate* and if we want to integrate further we have to call *Iterate* again but with a larger value of time.

For the sake of brevity I gonna create now a GraphicsGrid iterating through the *StateData* structure:

    GraphicsGrid[Partition[sols = Table[
    NDSolve`Iterate[stateD, t];
    Plot3D[
     Evaluate[
      Abs@u[t, x, y] /.
       NDSolve`ProcessSolutions[stateD, "Forward"]],
     {x, -5, 5}, {y, -5, 5}, PlotRange -> {0, 1}],
    {t, 0, 2, .1}], 2]]


![enter image description here][5]

And here the canonical animation:

![enter image description here][6]

This solution may not be the "straightforward" one, but if you are short on ressources (memory) this is the way to go.


According to *Weierstrass*, the ultimative goal is always the representation of a function.

I guess I fulfilled his dictum now...

I hope this helps.


  [1]: http://i.stack.imgur.com/U3KDt.gif
  [2]: http://i.stack.imgur.com/jfuye.png
  [3]: http://i.stack.imgur.com/WcWVL.png
  [4]: http://reference.wolfram.com/mathematica/tutorial/NDSolveStateData.html#93858351%E2%80%8C%E2%80%8B
  [5]: http://i.stack.imgur.com/BwREb.png
  [6]: http://i.stack.imgur.com/qzxaL.gif
