Dear Gyoto developers,

A few notes on maintenance:

Before each commit:
===================

* If you touched `configure.ac` or any `Makefile.am` file, run `autoreconf`
  and commit the regenerated files (`configure`, `Makefile*.in`)

* Run `make check-all`. If possible, build the Python bindings first, so
  that the test suite also runs the detailed Python regression tests.

* Check that every file you are committing carries correct Copyright
  information, including the years. If in doubt, run
  `git blame <filename>` to see who modified the file and when.

* Add a line to the `NEWS` file, in the section for the upcoming (not yet
  released) version.

* If this commit changes the user manual, update the date on the title
  page.

To add a new class:
===================

* Create `lib/NewClass.C` and `include/GyotoNewClass.h`, e.g. by copying an
  existing file with a similar architecture and adapting it.
  See "Coding rules and caveats" below.

* Update `lib/Makefile.am` by adding your NewClass where it fits
  (typically in the list `libgyoto_stdplug_la_SOURCES`)

* Update `lib/StdPlug.C` (or the corresponding file of the plug-in in
  which your class should appear): add the two lines that register your
  class, copying what the neighbouring classes do.

* Update `python/include/gyoto_std.i` and `python/include/gyoto_std_swig.h`
  to declare your class, following the pattern used by similar classes.

* Run `autoreconf`

* Run `configure` and `make`

* Test (run `make check-all!`), debug, add an XML example file in
  `doc/examples/`, and consider adding a few lines to `python/example.py`
  or a test case to `python/tests/test_std.py`

* Push:
  - `lib/Makefile.am`
  - `lib/Makefile.in`
  - your `NewClass.C`
  - your `GyotoNewClass.h`
  - `lib/StdPlug.C`
  - `python/include/gyoto_std.i`
  - `python/include/gyoto_std_swig.h`
  - `doc/examples/example-newclass.xml`

To make a new release:
======================

* Check with the team whether anyone has plans for this release.

* Update the date in the user manual. Use the date of the last change
  to the manual, as shown by `git log doc/user_guide`, and edit
  `doc/user_guide/GyotoManual.tex` accordingly.

* Choose a version number. The version is MAJOR.MINOR.MICRO. Bump MAJOR
  when significant new features are added; bump MINOR for continued
  development of existing features; bump MICRO for bugfix or packaging
  releases.

* Edit `configure.ac`
  + update the version info in `AC_INIT`
  + update `VERSINFO` according to
    <http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html>
    Note that this `VERSINFO` is the libtool library version and has
    nothing to do with the Gyoto version.
  + run `autoreconf`.

* Edit `NEWS`
  + add a new entry (or edit the one that has not been released yet)
  + mark the release as "ABI n" if the ABI changed, "BUG" if it is a
    backwards-compatible bugfix release, "PKG" if only the packaging
    changed.
  + briefly document what changed in Gyoto: new features, significant
    bug fixes...

* Regenerate ChangeLog:
  + run `git log --no-merges > ChangeLog`
  + add an item on top corresponding to the log message you will use
    for the *next* commit (the one that performs the release). Use
    `unknown` as the commit identifier and the output of
    `LANG=C date --rfc-822` as the date.

* Clean rebuild and full test: `make distclean`, `./configure`, `make`,
  `make check`.

* Commit the changes.

* If releasing a branch other than stable: merge this branch into
  stable, make sure the two branches merge cleanly in both directions
  and are identical, then merge stable into master.

* Tag the release using the version number as the tag name, with a
  signed tag: `git tag -s 1.0.0` (replace 1.0.0 with the actual
  version; this requires a GPG key).

* Take a coffee break.

Coding rules and caveats:
=========================

Strings and numbers:
--------------------
Never use `std::atof`, `std::stod` or similar generic functions to convert
between strings and doubles: they obey the user's locale (on a French
system, "3.14" may fail to parse). Use `Gyoto::stringToDouble` and
`Gyoto::doubleToString` instead: they always use "." as the decimal
point and handle special values such as `DBL_MAX`.

The same applies to output, especially XML files, which must always use
"." as the decimal point whatever the user's locale is. If using
C-style `sprintf()`, do:

```cpp
#include <clocale>
std::string loc(setlocale(LC_NUMERIC, NULL));
setlocale(LC_NUMERIC, "C");
sprintf(txt, fmt, val);
setlocale(LC_NUMERIC, loc.c_str());
```

If using C++ streams, do:

```cpp
#include <locale>
std::ostringstream ss;
ss.imbue(std::locale("C"));
ss << val;
```

Overriding methods in a derived class (recipe):
------------------------------------------------
Gyoto is built on inheritance: a base class such as
`Gyoto::Metric::Generic` (in `GyotoMetric.h`) declares methods that
derived classes such as `Gyoto::Metric::KerrBL` redefine. Three C++
terms are involved; you do not need to master them, just follow the
recipe below:

* "virtual" (written in the base class) marks a method that derived
  classes may redefine.
* "override" (written in your derived class) tells the compiler "this
  method replaces a base-class method". If the argument list does not
  match exactly, you get a compile-time error instead of a silent bug.
* "overloaded" means the base class provides several methods with the
  same name but different arguments, e.g. Metric Generic::gmunu exists
  as both `gmunu(double const x[4], int mu, int nu)` and `gmunu(double
  g[4][4], double const x[4])`. If you redefine only one of them, C++
  hides the others: calls to the hidden versions will fail.

The recipe, to be applied in your class declaration:

```cpp
class KerrBL: public Generic {
  ...
public:
  // Rule 1: un-hide the overloads we don't redefine
  using Generic::gmunu;
  // Rule 2: mark the one we do redefine with 'override'
  double gmunu(double const x[4], int mu, int nu) const override;
};
```

Rule 1: always add the `using BaseClass::method;` line for every
overloaded base-class method you redefine.
Rule 2: always end your redefined method with the `override` keyword.

Gyoto builds with `-Woverloaded-virtual -Werror=overloaded-virtual` in
`CXXFLAGS`, so breaking Rule 1 is a hard compile-time error: the
compiler tells you exactly which method is hidden and where to add the
"using" line.

Following Rule 2 ensures that a typo in the parameter list of a method
meant to override a base-class method — or a change in the base class
not yet reflected in the derived class — produces an explicit
compile-time error instead of a hard-to-find bug.

Documenting your class:
-----------------------
* Document the class and its non-inherited methods in the header file
  (`include/GyotoNewClass.h`) with Doxygen comments. This text is
  published on <https://gyoto.obspm.fr/> and in the Python help:
  `help(gyoto.std.NewClass)` at the Python prompt.

* Document each custom property with the `GYOTO_PROPERTY_*` macros in the
  `.C` file. This description shows up in Python (`myinstance.help()`)
  and when hovering widgets in graphical user interfaces.

Regards, Thibaut.
