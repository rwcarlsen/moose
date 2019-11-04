
system=$1

perl -0777 -i -pe 's/template *<>\s*InputParameters\s*validParams<([A-Za-z0-9_]*)>\(\)/defineLegacyParams($1);\n\nInputParameters\n$1::validParams()/gs;' -pe 's/(InputParameters|auto) (\w+) = validParams<([a-zA-Z0-9_]*)>\(\)/$1 $2 = $3::validParams\(\)/gs;' framework/src/$system/*.C

perl -0777 -i -pe 's/(class.*: public .*{\s*public:\s*)/$1static InputParameters validParams();\n\n  /gs;' framework/include/$system/*.h
