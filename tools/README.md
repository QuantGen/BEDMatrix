This folder contains supporting files for Autoconf.

The macros `m4/ax_boost_base.m4` and `m4/ax_boost_iostreams.m4` are used to
detect Boost in the configure stage and were downloaded from the [GNU Autoconf
Archive](http://www.gnu.org/software/autoconf-archive/).

`config.guess`, `config.sub`, and `install-sh` are dependencies of
`m4/ax_boost_base.m4`. `config.guess` and `config.sub` were obtained using
instructions from the [GNU gettext
manual](https://www.gnu.org/software/gettext/manual/html_node/config_002eguess.html),
and `install-sh` was copied from a local copy of Automake 1.15.
