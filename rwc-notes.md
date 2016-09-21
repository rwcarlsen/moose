
Thoughts/questions
=====================

* Kernel API should pass in quadrature point value (_qp) to computeQpResidual,
  computeQpJacobian, etc.  It should also pass in _i and possible other
  values.  It might be worth gathering them into a struct.

* Overloading operators is confusing (IMO).  I think calling member functions
  explicitly is more clear/better.

* Auto-formatting issue - I think that it would be best to adapt code style
  guide to something that can be clang-formatted.  Since style is arbitrary
  (given most reasonable styles), the time saved will quickly pay itself off
  beyond any downsides.  Hand-rolling tools to accomodate current style
  details will be a waste of time.  Style checking (i.e. not reformatting)
  tools add more complexity to process and still require lots of manual code
  inspection and modification.

* I don't like the code duplication required by users in kernels - i.e.:
    
```
template<>
InputParameters validParams<ExampleConvection>()
{
  InputParameters params = validParams<Kernel>();
  // the line here:
  params.addRequiredParam<RealVectorValue>("velocity", "Velocity Vector");
  return params;
}

ExampleConvection::ExampleConvection(const InputParameters & parameters) :
  // You must call the constructor of the base class first
  Kernel(parameters),
   // is quite duplicated here:
   _velocity(getParam<RealVectorValue>("velocity"))
{}
```

* Was the idea of having a new app-generation tool considered?  Instead of a
  repository that you fork, we could have a tool where the user passes a name
  for their kernel/app and we populate a directory full of build files,
  header/implementation files, etc. automatically with everything renamed
  appropriately.  Downside is don't get to track usage via github forking
  record and also makes it possible for users to experiment without using
  version control (could be both good and bad).

* Boost-serialization for checkpointing? - it automagically supports recursive
  stl containers, etc. into a nice human readable XML format.  Checkpointing
  usually doesn't need to be fast - it is really only needed occasionally
  during runs.

* time-interval-based checkpointing (would require threading).

* Why can only scalar material properties have defaults?

* Document output file format(s) and how to manipulate them - are there cli
  utils?

* Support provenance with storing meta-data about a run - e.g. moose version,
  optional deps installed, deps versions, input file, run config (e.g. mpi,
  etgc.), etc.  Maybe keep it all together too in some sort of archive/db.

* Already use yaml for some stuff - has it been considered as main moose input
  file format?

* Why does MooseArray exist? - why not just std::vector?

* Why type 'Real' - and where is it defined?

* Why do we expose ``_test`` and ``_grad_test`` vectors to user - why not just
  set them directly as the current functions (i.e. pre-index with ``_i`` and
  ``_j``)?

* Put all moose code under Moose namespace

