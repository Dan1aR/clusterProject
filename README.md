# clusterProject
Проект по кластеризации, Андреев Данил ИУ6-23Б, МГТУ им. Н.Э.Баумана. 
Алгоритмы: K-средних, Иерархический, Волновой

Описание:

pointsCreator.cpp - программа, создающая кластера по следующему шаблону:
(x, y, номер кластера к котрому принадлежит точка)
Данные сохраняются в input.txt, он будет приложен для тестирования.
Номер кластера не используется в main.cpp, он нужен для визуализации input.txt

main.cpp - программа, осуществляющая кластеризацию
Она пока не содержит интерфейса, только три алгоритма кластеризации 
(в main() вызывается один, остальные два из них закомментированы, чтобы не мешаться)
Глобальные переменные, которые должен вводить пользователь, определены сразу после команд препроцессора.
K = 3 - количество кластеров (для k-means)
INTERVAL = 50 - интервал сетки (для waveCluster)
FILENAME = "input.txt" - файл, из которого считываются точки

результат сохраняется в output.txt по такому же шаблону

# gnuplot
Рекомендую следующую команду для визуализации:
plot "< awk '{if($3 == \"0\") print}' output.txt" u 1:2 t "red" w p pt 2, "< awk '{if($3 == \"1\") print}' output.txt" u 1:2 t "green" w p pt 2, "< awk '{if($3 == \"2\") print}' output.txt" u 1:2 t "blue" w p pt 2, "< awk '{if($3 == \"3\") print}' output.txt" u 1:2 t "pink" w p pt 2
