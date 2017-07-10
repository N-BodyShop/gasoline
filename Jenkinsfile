pipeline {

  agent any

  stages {

    stage('Checkout') {
      steps {
        checkout scm
        sh 'cp Makefile.default Makefile'
      }
    }

    stage('Build') {
      steps {

        echo 'building gasoline with pthreads...'
        sh 'make EXE="gasoline_pthread" pthread'
        sh 'make clean'

        echo 'building gasoline with mpi...'
        sh 'module load mpi'
        sh 'make EXE="gasoline_mpi" mpi'
        sh 'make clean'

        sh 'rm Makefile'
      }
    }
  }

  post {
    success {
      echo 'build success'
    }
    failure {
      echo 'build failure'
    }
  }
}
