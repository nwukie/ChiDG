node {
    stage('Checkout'){
        git url: 'https://github.com/nwukie/ChiDG.git', branch: 'dev'
    }

    sh 'mkdir -p build'
    dir('build'){    

        
        stage('Configure'){
                sh('cmake -DCMAKE_Fortran_COMPILER=gfortran ..')
            }
    
        stage('Build Core'){
            sh 'make -j 2'
        }
    
        stage('Build Test'){
            sh 'make check -j 2'
        }
    
        stage('Test - pFUnit'){
            sh 'make test'
            junit testResults: 'bin/*.xml' 
        }
        
    }
}
