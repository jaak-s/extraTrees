Installation for development
============================

1. clone git repository and import it as Java project in Eclipse
2. run ant (in Eclipse) with the included *pom.xml* and make *dist*
  * for configuring ant, add Environment variable (External Tools Conf -> Environment)
     * Key: JAVA6_BOOT_CLASSPATH
     * Val: /usr/lib/jvm/java-6-openjdk/jre/lib/rt.jar

Ant creates a folder "~/Documents/temp/extraTrees" which can be installed into R:
```bash
cd ~/Documents/temp
R CMD INSTALL extraTrees
## to build tar.gz of the package:
R CMD build extraTrees  
```

The previous steps assume you have **rJava** installed.

## Windows
In windows use following command to install the package
```R
install.packages("C:\\Users\\...\\extraTrees", type="source", repos=NULL)
```

## Running included tests
After the package is installed you can run the R tests.
```R
library(testthat)
test_package("extraTrees")
```
 
