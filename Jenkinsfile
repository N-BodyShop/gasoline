pipeline {

  agent any

  stages {

    stage('Checkout') {
      steps {
        checkout scm
        sh 'cp Makefile.default Makefile'
      }
    }

    stage('Build pthread') {
      steps {
        echo 'building gasoline with pthreads...'
        sh 'make EXE="gasoline_pthread" pthread'
        sh 'make clean'
      }
    }

    stage('Build mpi') {
      steps {
        echo 'building gasoline with mpi...'
        sh ". /etc/profile.d/modules.sh \n" +
           "module load mpi \n" +
           "make EXE=\"gasoline_mpi\" mpi"
        sh 'make clean'
      }
    }

  }

  post {
    always {
      deleteDir()
    }
    success {
      echo 'build success'
    }
    failure {
      echo 'build failure'
    }
  }

  options {
    timeout(time: 60, unit: 'MINUTES')
  }
}
