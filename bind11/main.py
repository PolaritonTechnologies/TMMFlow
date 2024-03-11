import example

# Call the function
result = example.add(2, 3)
print(result)  # Outputs: 5

#  c++ -O3 -Wall -shared -std=c++17 -fPIC -IC:\Users\GatherLab-Julian\AppData\Local\Programs\Python\Python39\lib\site-packages\pybind11\include -IC:\Users\GatherLab-Julian\AppData\Local\Programs\Python\Python39\Include example.cpp -o example.pyd
#  g++ -O3 -Wall -shared -std=c++17 -fPIC -IC:\Users\GatherLab-Julian\AppData\Local\Programs\Python\Python312\Lib\site-packages\pybind11\include -IC:\Users\GatherLab-Julian\AppData\Local\Programs\Python\Python312\Include -LC:\Users\GatherLab-Julian\AppData\Local\Programs\Python\Python312\libs -lpython312 example.cpp -o example.pyd