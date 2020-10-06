        program test

        double precision alpha

5       write(6,10,advance='no')
10      format(' Enter alpha: ')
        read(5,*) alpha
        write(6,*) cos(alpha),sin(alpha)
        go to 5

        end 
