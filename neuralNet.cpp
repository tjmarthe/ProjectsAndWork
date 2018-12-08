//neuralNet Implemented by Trevor Marthe

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>


struct Connection
{
  double weight;
  double weightDiff;
};

class Neuron;

typedef std::vector<Neuron> Layer;

class Neuron {

public:

  double transferFunc(double num){
    return num / (1 + abs(num));
  }

  double transferFuncDerv(double num){
    return  num * (1-num);
  }

  Neuron(int numOut, int nerIndex_){
    for(int i = 0; i < numOut; ++i){
      outputWeights.push_back(Connection());
      outputWeights.back().weight = rand() / double(RAND_MAX);
    }
    nerIndex = nerIndex_;
  }
  void setOut(double num){ outputVal = num;}
  double getOut(){return outputVal;}
  void feedForward(Layer prevlayer){
    //std::cout << "feed ner 1" << std::endl;
    double result;
    for(int ner = 0; ner < prevlayer.size(); ++ner){
      //      std::cout << "feed ner one" << std::endl;
      result += prevlayer[ner].getOut() * prevlayer[ner].outputWeights[nerIndex].weight;
      //std::cout << "feed ner" << std::endl;
    }
    outputVal = transferFunc(result);
  }


  double sumDOW(Layer nextlayer){

    double result;

      for (int ner = 0; ner < nextlayer.size() -1; ++ner){
        result += outputWeights[ner].weight + nextlayer[ner].grad;

      }

    return result;
  }

void calcOutGrad(double num){

  double diff = num - outputVal;

  grad = diff * transferFuncDerv(outputVal);
}

void calcHidGrad(Layer nextlayer){

  double result = sumDOW(nextlayer);
  grad = result * transferFuncDerv(outputVal);

}


void updateInWeights(Layer prevlayer){

  for(int ner = 0; ner < prevlayer.size(); ++ner){
    Neuron neuron = prevlayer[ner];
    double oldDiffWeight = neuron.outputWeights[nerIndex].weightDiff;

    double newDiffWeight = learningRate * neuron.getOut() * grad + alpha * oldDiffWeight;

    neuron.outputWeights[nerIndex].weightDiff = newDiffWeight;
    neuron.outputWeights[nerIndex].weight += newDiffWeight;
  }

}

private:

double grad; //gradient
int nerIndex;
double outputVal;
std::vector<Connection> outputWeights;
double learningRate = .20; //changable
double alpha = .3;
};


class Net{
public:
  Net(std::vector<unsigned> construct){

    unsigned numofLayers = construct.size();
    for(unsigned i = 0; i < numofLayers; ++i){

      layers.push_back(Layer());
      unsigned numofOuts = i == construct.size() - 1 ? 0 : construct[i + 1];

      for(unsigned ner = 0; ner <= construct[i]; ++ner){
        layers.back().push_back(Neuron(numofOuts, ner));
      }
    }

    layers.back().back().setOut(1.0); //setting bias of the last layer
    
    

  }

  void feedForward(std::vector<double> inputVals){
    for (int i = 0; i < inputVals.size(); ++i){
      //      std::cout << "layer" << std::endl;
      layers[0][i].setOut(inputVals[i]);
    }

    //forward prop
    for(int i = 1; i < layers.size(); ++i){
      //std::cout << "feeding" << std::endl;
      Layer &prevLayer = layers[i -1];//for the neuron to know the prevlayer
      //std::cout << "prev layer" << std::endl;
      for(int ner = 0; ner < layers[i].size()-1; ++ner){//minus 1 cause bias neuron
        //std::cout << "feed forward" << std::endl;
        layers[i][ner].feedForward(prevLayer);
        //std::cout << "feed men" << std::endl;
      }
    }


  }


  void backProp(std::vector<double> targets){
    Layer outLayer = layers.back();
    error = 0.0;

    for (int ner = 0; ner < outLayer.size() - 1; ++ner){
      double diff = targets[ner] - outLayer[ner].getOut();
      error += diff * diff;
    }

    error /= outLayer.size() - 1;
    error = sqrt(error);

    std::cout << error << std::endl;


    recentAErr = (recentAErr* recentASmoothErr + error)/ (recentAErr+1.0);

    for (int ner = 0; ner < outLayer.size() -1; ++ner){
      outLayer[ner].calcOutGrad(targets[ner]);
    }

    for (int i = layers.size()-1; i > layers.size() - 0; --i)
      {
        Layer hiddenLayer = layers[i];
        Layer nextLayer = layers[i+1];

        for(int ner = 0; ner < hiddenLayer.size(); ++ner)
          hiddenLayer[ner].calcHidGrad(nextLayer);
      }


    for (int ner = 0; ner < outLayer.size() -1; ++ner){
      outLayer[ner].calcOutGrad(targets[ner]);
    }


    for (int i = layers.size()-1; i > layers.size() - 0; --i)
      {
        Layer Layer1 = layers[i];
        Layer prevLayer = layers[i+1];

        for(int ner = 0; ner < Layer1.size()-1; ++ner)
          Layer1[ner].updateInWeights(prevLayer);
      }



  }

  void getResult(std::vector<double> resultVals){
    resultVals.clear();
    //std::cout << "clear" << std::endl;
    for(int i = 0; i < layers.size(); ++i){
      //      std::cout << layers.back()[i].getOut() << " output actual" << std::endl;
      resultVals.push_back(layers.back()[i].getOut());
    }
    std::cout<< "OUTPUT: ";
    for(int x = 0; x < resultVals.size(); ++x){
      if(resultVals[x] > .750)
        std::cout << 1 << " ";
      else
        std::cout << 0 << " ";
    }
    std::cout << std::endl;
}

private:

  std::vector<Layer> layers;
  double error;
  double recentAErr;
  double recentASmoothErr;

};



int main(){

  srand(time(0));

  std::vector<unsigned> construct;

  construct.push_back(4);
  construct.push_back(5);
  construct.push_back(4);

  Net myNet(construct);


  std::vector<double> input;
  std::vector<double> output;
  std::vector<double> target;
  double input1, input2, input3, input4;
  //double target1, target2, target3;

  /*
  input is a simulation of  button pressing,
  a combination of 4 buttons is pressed, based
  on that combination, the out put of 3 lights
  wll change.
  ie: if buttons 1 & 2(input 1 and 2 aree both 1) are pressed, then light 1(output 1 is 1) will ligh\
t up.
      if buttons 3 & 4(input 3 and 4 are both 1) are pressed, then light 3(output 3 is 1) will light\
 up.
      All buttons = All outputs
      2&3 = output 2 is 1
  */

  for(int i = 0;  i < 100; ++i){
    input.clear();
    output.clear();
    target.clear();
    input1 = double(rand()) / double(RAND_MAX);;
    //    std::cout << input1 << std::endl;
    if(input1 < .49)
      input.push_back(0.0);
    else
      input.push_back(1.0);
    input2 = double(rand()) / double(RAND_MAX);;
    if(input2 < .49)
      input.push_back(0.0);
    else
      input.push_back(1.0);
    input3 = double(rand()) / double(RAND_MAX);;
    if(input3 < .49)
      input.push_back(0.0);
    else
      input.push_back(1.0);
    input4 = double(rand()) / double(RAND_MAX);;
    if(input4 < .49)
      input.push_back(0.0);
    else
      input.push_back(1.0);
    std::cout << "INPUT: " ;
    for(int j = 0; j < input.size(); ++j)
      std::cout << input[j] <<  " ";
    std::cout << std::endl;
    for(int IT = 0; IT < input.size(); ++IT){
      //std::cout << "target" << std::endl;
      if(IT == input.size()-1)
        break;
      if(input[IT] == 0.0 && input[IT + 1] == 0.0)
        target.push_back(0.0);
      if(input[IT] == 1.0 && input[IT + 1] == 1.0)
        target.push_back(1.0);
      else
        target.push_back(0.0);
    }

    //std::cout<< "input" << std::endl;
    myNet.feedForward(input);
    //std::cout << "get out" << std::endl;
    
     myNet.getResult(output);
    myNet.backProp(target);


  }
}




