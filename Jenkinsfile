node {
    stage('Checkout'){
        checkout([$class: 'GitSCM', branches: [[name: 'dev']], doGenerateSubmoduleConfigurations: false, extensions: [], submoduleCfg: [], userRemoteConfigs: [[credentialsId: '31bb08f6-d8c2-47f0-a217-af319cc191f9', url: 'https://github.com/nwukie/ChiDG.git']]])
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
