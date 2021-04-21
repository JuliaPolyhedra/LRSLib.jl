using BinDeps
using Libdl

@BinDeps.setup

lrslib_commit = "5b8ab79d3b3a40f128083da3935523cb254d3460"
#lrsname = "lrslib-061"
lrsname = "lrslib-$lrslib_commit"
lrsnashname = "lrsnashlib-$lrslib_commit"
lrslibname = "liblrs"

# julia installs libgmp10 but not libgmp-dev since it
# does not have to compile C program with GMP,
# e.g. it does not need the headers.
# FIXME BinDeps doesn't work with header, it only works for libraries because it checks for the .so
#       It sees libgmp.so (installed by libgmp10 as a julia dependency) and thinks that it's ok but
#       it does not have the headers
#libgmpdev = library_dependency("libgmp-dev", aliases=["libgmp"])
liblrs = library_dependency("liblrs", aliases=[lrsname, "liblrsgmp"])#, depends=[libgmpdev])
liblrsnash = library_dependency("liblrsnash", aliases=[lrsnashname, "liblrsnashgmp"])

official_repo = "http://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/$lrsname.tar.gz"
forked_repo = "https://github.com/JuliaPolyhedra/lrslib/archive/$lrslib_commit.zip"

#GMP
@static if Sys.islinux()
    const has_apt = try success(`apt-get -v`) catch e false end
    const has_yum = try success(`yum --version`) catch e false end
    const has_pacman = try success(`pacman -Qq`) catch e false end
    if has_apt || has_yum
        if has_apt
            pkgname = "libgmp-dev"
            pkgman = "apt-get"
        else
            pkgname = "libgmp-devel or gmp-devel"
            pkgman = "yum"
        end

        println("Warning: The compilation of LRS requires the header gmp.h provided by the package $pkgname.")
        println("If the compilation fails, please install it as follows:")
        println("\$ sudo $pkgman install $pkgname")
    end
end

#LRS
provides(Sources,
         Dict(URI(forked_repo) => liblrs), unpacked_dir="$lrsname")

lrssrcdir = joinpath(BinDeps.srcdir(liblrs), lrsname)
lrsprefixdir = joinpath(BinDeps.usrdir(liblrs))
lrslibdir = joinpath(lrsprefixdir, "lib")

targetdirs = AbstractString["$lrslibname.$(Libdl.dlext)"]
patchdir = BinDeps.depsdir(liblrs)

@static if Sys.isapple()
    using Homebrew
    Homebrew.add("gmp")
    homebrew_includedir = joinpath(Homebrew.brew_prefix, "include")
    homebrew_libdir = joinpath(Homebrew.brew_prefix, "lib")
    provides(BuildProcess,
             (@build_steps begin
              GetSources(liblrs)
              CreateDirectory(lrsprefixdir)
              CreateDirectory(lrslibdir)
              @build_steps begin
              ChangeDirectory(lrssrcdir)
              FileRule(joinpath(lrslibdir,"$lrslibname.$(Libdl.dlext)"),@build_steps begin
                       pipeline(`cat $patchdir/makefile.osx.patch`, `patch`)
                       `make $lrslibname.$(Libdl.dlext).0 SONAME=$lrslibdir/$lrslibname.$(Libdl.dlext) SHLIB=$lrslibname.$(Libdl.dlext).0 SHLINK=$lrslibname.$(Libdl.dlext) INCLUDEDIR=$homebrew_includedir LIBDIR=$homebrew_libdir`
                       `cp $lrslibname.$(Libdl.dlext).0 $lrslibdir/$lrslibname.$(Libdl.dlext)`
                       end)
              end
             end),liblrs)
    provides(BuildProcess,
             (@build_steps begin
              ChangeDirectory(lrssrcdir)
              FileRule(joinpath(lrslibdir,"liblrsnash.$(Libdl.dlext)"),@build_steps begin
                       `cp $patchdir/makefile.liblrsnash $lrssrcdir`
                       `make -f makefile.liblrsnash SHLIB=liblrsnash.$(Libdl.dlext) LIBLRSDIR=$lrslibdir INCLUDEDIR=$homebrew_includedir LIBDIR=$homebrew_libdir`
                       `cp liblrsnash.$(Libdl.dlext) $lrslibdir/liblrsnash.$(Libdl.dlext)`
                       end)
             end),liblrsnash)
else
    provides(BuildProcess,
             (@build_steps begin
              GetSources(liblrs)
              CreateDirectory(lrsprefixdir)
              CreateDirectory(lrslibdir)
              @build_steps begin
              ChangeDirectory(lrssrcdir)
              FileRule(joinpath(lrslibdir,"$lrslibname.$(Libdl.dlext)"),@build_steps begin
                       `make $lrslibname.$(Libdl.dlext).0 SONAME=$lrslibname.$(Libdl.dlext) SHLIB=$lrslibname.$(Libdl.dlext).0`
                       `cp $lrslibname.$(Libdl.dlext).0 $lrslibdir/$lrslibname.$(Libdl.dlext)`
                       end)
              end
             end),liblrs)
    provides(BuildProcess,
             (@build_steps begin
              ChangeDirectory(lrssrcdir)
              FileRule(joinpath(lrslibdir,"liblrsnash.$(Libdl.dlext)"),@build_steps begin
                       `cp $patchdir/makefile.liblrsnash $lrssrcdir`
                       `make -f makefile.liblrsnash SHLIB=liblrsnash.$(Libdl.dlext) LIBLRSDIR=$lrslibdir`
                       `cp liblrsnash.$(Libdl.dlext) $lrslibdir/liblrsnash.$(Libdl.dlext)`
                       end)
             end),liblrsnash)
end

BinDeps.@install Dict([(:liblrs, :liblrs),
                       (:liblrsnash, :liblrsnash)])
