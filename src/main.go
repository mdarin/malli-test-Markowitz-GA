//
// Mally lib test and usage example
// This program creates optimal portfolio(allocate assets) 
// using genetic algorithm to find the most attractive 
// Risk/Return ratio by Markowitz model
//
package main

import(
	// common purpose
	"fmt"
	"os"
	_ "log"
	"time"

	// mally platfom
	ma "github.com/mdarin/mally"

	// utils
	"sort"
	"sync"
	"math"
)






// Интерфейс к портфелю
type PortfolioSInterface interface {
	// интерфейс для методов оптимизации
	// функция g(x) оценки минимума в методах одномерной оптимизации 
	// в методе наискорейшего спуска(вообще в многомерной оптимизации) это функция g(x[k]) 
	// применяемая для нахождения шага поиска lambda[k] как минимума функции g(x[k])
	// на отрезке a,b 
	// lambda здесь это вычисленная точка по методу поиска минимума, имя надо бы выбрать другое... 
	GetG() func(data ma.ODO, lambda float64) float64
	// фукнция миниму которой мы ищем в методах многомерной оптимизации
	Function() float64
	SetFunction(f func(float64) float64)
	// вектор аргументов функции
	GetX() []float64
	SetX(args []float64)
	// значение фукнции
	SetY(float64)
	GetY() float64
	// градиет функции - вектор частных производных
	Grad() []float64
	// отрезок для поиска минимума методом одной оптимизации
	//a,b float64
	GetInterval() (float64, float64)
	// максимальное количесво приближений 
	//n int - для одномерной оптимизации
	GetMaxIterODO() int
	//m int - многомерной оптимизации
	GetMaxIterMDO() int
	// требуемая точность
	//epsilon1 float64 для метода одномерной оптимизации
	GetEpsilon1() float64
	//epsilon2 float64 для метода многомерной оптимизации
	GetEpsilon2() float64
	// Аргумент минимизации (для наискорейшего спуска argmin по lambda)
	// phi(x) метод одномерной оптимизаци
	GetF() (func(data ma.ODO) float64)

	// методы опосредованного достпу к элементам описывающим портфель
	GetReturn() float64
	SetReturn(ret float64)
	GetRisk() float64
	SetRisk(risk float64)
	GetWeights() *ma.Matrix
	SetWeights(w *ma.Matrix)
	GetAssets() *ma.Matrix
	SetAssets(a *ma.Matrix)
	GetTickers() *ma.Matrix
	GetRiskReturnRatio() float64

	GetTempMeanReturns() *ma.Matrix
	GetTempCovariances() *ma.Matrix
	GetTempA() *ma.Matrix
	GetTempAinv() *ma.Matrix
	GetTempB() *ma.Matrix
	GetTempX() *ma.Matrix
}

//Вроде такой структура портфеля должна быть
type PortfolioS struct {
	tickers *ma.Matrix // - вектор индексов основной матрицы указывает на столбцы с данными
	assets *ma.Matrix // - матрица(таблица) инструментов
	weights *ma.Matrix // - вектор весов
	risk float64 //- дисперсия портфеля
	ret float64 //- доходность
	// минимальный анализ:
	sigma float64 // - СКО портфеля (корень из дисперсии)
	retMin float64 //= Return - 3Sigma
	retMax float64 //= Return + 3Sigma 
	// служебные поля
	nrows int
	ncols int
	// временные таблицы(матрицы)
	rx *ma.Matrix
	cov *ma.Matrix
	a *ma.Matrix
	a_inv *ma.Matrix
	b *ma.Matrix
	x *ma.Matrix
}
// constructor
func New(nrows, ncols int) (*PortfolioS,error) {
	var err error = nil
	var p *PortfolioS = new(PortfolioS)
	p.tickers,err = ma.New(ncols,1) // vector col
	p.assets,err = ma.New(nrows,ncols)
	p.weights,err = ma.New(ncols,1) // vector col
	// создаём временные таблицы(матрицы)
	p.rx,err = ma.New(ncols, 1)
	p.cov,err = ma.New(ncols, ncols)
	p.a,err = ma.New(p.cov.Rows()+2, p.cov.Cols()+2)
	p.a_inv,err = ma.New(p.a.Rows(), p.a.Cols())
	p.b,err = ma.New(p.a.Rows(), 1)
	p.x,err = ma.New(p.a.Rows(), 1)
	p.nrows = nrows
	p.ncols = ncols
	return p,err
}
func (p *PortfolioS) Rows() int { return p.nrows }
func (p *PortfolioS) Cols() int { return p.ncols }

// portfolio operators
func (p *PortfolioS) Make(a *ma.Matrix) error {
	//
	// Создаём и инициализируем портфель
	//
	var err error
	// сформировать портфель
	// генерируем вектор индексов по желаемому кочиству бумаг в портфеле
	// так как распределение равномерное, беспокоиться не надо, все акции 
	// поучаствуют в отборе практически равновероятно
	fmt.Printf("portfolio: \n")
	// выбрать индкес ячейки в портфеле - номер колонки 0,1,2,3... 
	// выбрать индеск тикера из массива - номер колонки
	// построчно скопировать данные из массива в портфель
	// пока ячейки в портфеле не закончатся
	fmt.Printf("  tickers: ")
	for j := 0; j < p.GetAssets().Cols(); j++ {
		// собрать портфель
		//TODO:
		ticker := ma.RndBetweenU(0, a.Cols())
		//FIXME: DEBUG!!!!!!
		//ticker := j // DEBUG

		fmt.Printf("[%d]", ticker)
		p.GetTickers().Setij(j, 0, float64(ticker))
		for i := 0; i < p.GetAssets().Rows(); i++ {
			value,_ := a.Getij(i, ticker)
			p.GetAssets().Setij(i, j, value)
		}
	}
	fmt.Printf("\n")
	//ma.Mprintf("#Assets:", " %.3f", p.GetAssets())

	return err
}

func (p *PortfolioS) Formalize() {
	var workers int
	var wg sync.WaitGroup

	//
	// Произвести формулировку задачи оптимизации инвестиционного портфеля
	// ценных бумаг по Маркивицу
	//

	// получаем временные таблицы(матрицы)
	Assets := p.GetAssets()
	Rx := p.GetTempMeanReturns()
	COV := p.GetTempCovariances()
	A := p.GetTempA()
	Ainv := p.GetTempAinv()



	// Вычисляем средние доходности инструментов
	workers = Assets.Cols()
	wg.Add(workers)
	for i := 0; i < workers; i++ {
		go func (a, b *ma.Matrix, cursor int, wg *sync.WaitGroup) {
			defer wg.Done()
			mean := ma.Mean(b, cursor, ma.COLVECTOR)
			a.Setij(cursor, 0, mean)
		}(Rx, Assets, i, &wg)
	} // eof for
	wg.Wait()

	// вычисляем ковариации для всех инструментов
	workers = Assets.Cols()
	wg.Add(workers)
	for i := 0; i < workers; i++ {
		go func (a, b *ma.Matrix, cursor int, wg *sync.WaitGroup) {
			defer wg.Done()
			for j := 0; j < b.Cols(); j++ { // тут ошибка было Rows() вместо Cols()!!!
				a.Setij(cursor, j, ma.Covariance(b, b , cursor, j, ma.COLVECTOR))
			}
		}(COV, Assets, i, &wg)
	} // eof for
	wg.Wait()

	// создаём матрицу системы уравнений по методу Лагранжа
	// craete matrix for inversing
	for i := 0; i < COV.Rows(); i++ {
		for j := 0; j < COV.Cols(); j++ {
			value,_ := COV.Getij(i, j)
			if i == j {
				A.Setij(i, j, 2 * value)
			} else {
				A.Setij(i, j, value)
			}
		}
	}
	for i := 0; i < A.Rows(); i++ {
		// 1
		A.Setij(i, A.Rows()-2, float64(1))
		A.Setij(A.Rows()-2, i, float64(1))
		// Ri
		value,_ := Rx.Getij(i, 0)
		A.Setij(i, A.Rows()-1, value)
		A.Setij(A.Rows()-1, i, value)
		// 0
		if i > A.Rows()-3 {
			A.Setij(i, A.Rows()-2, float64(0))
			A.Setij(A.Rows()-2, i, float64(0))
		}
	}

	// Получаем обратную матрицу для системы уравнений
	// при расчёте долей инструментов 
	// get A-1 matrix 
	epsilon := 0.000001 // праметр требуемой точности
	ma.Inverse(Ainv, A, epsilon)

	//ma.Mprintf("##Assets", " %f", Assets)
	//ma.Mprintf("\n#Rx", "%f", Rx)
	//ma.Mprintf("\n#COV", "%f ", COV)
	//ma.Mprintf("\n#A", "%f ", A)
	//ma.Mprintf("\n#Ainv", "%f ", Ainv)
}

// оператор присвоения :=
func (p *PortfolioS) Assign(value *PortfolioS) {
	//TODO: if value == nil return nil
	ma.Assign(p.tickers,value.tickers)
	ma.Assign(p.assets,value.assets)
	ma.Assign(p.weights,value.weights)
	ma.Assign(p.rx,value.rx)
	ma.Assign(p.cov,value.cov)
	ma.Assign(p.a,value.a)
	ma.Assign(p.a_inv,value.a_inv)
	ma.Assign(p.b,value.b)
	ma.Assign(p.x,value.x)
	p.risk = value.risk
	p.ret = value.ret
}



// требуется продумать операторы для портфеля(вроде готово)
// оператор включения в портфель актива из другого портфеля
// (производится замена одного актива другим)
func (p *PortfolioS) Replace(from *PortfolioS, index int) {
	// получить значение тикера ячейки с индексом index источника
	ticker,_ := from.GetTickers().Getij(index, 0)
	// устаночить знаение тикера в ячейку с индексом index приёмника
	p.GetTickers().Setij(index, 0, ticker)
	// скопировать значения временного ряда
	for i := 0; i < p.GetAssets().Rows(); i++ {
		value,_ := from.GetAssets().Getij(i, index)
		p.GetAssets().Setij(i, index, value)
	}
}
// оператор включения актива из доступного множества активов
func (p *PortfolioS) Take(assets *ma.Matrix, ticker, index int) {
	// устаночить знаение тикера в ячейки с индексом index приёмника
	p.GetTickers().Setij(index, 0, float64(ticker))
	// скопировать значения временного ряда из основной таблицы в портфель
	for i := 0; i < p.GetAssets().Rows(); i++ {
		value,_ := assets.Getij(i, ticker)
		p.GetAssets().Setij(i, index, value)
	}
}


// interface
func (p *PortfolioS) GetReturn() float64 { return p.ret }
func (p *PortfolioS) SetReturn(ret float64) { p.ret = ret }
func (p *PortfolioS) GetRisk() float64 { return p.risk }
func (p *PortfolioS) SetRisk(risk float64) { p.risk = risk }
func (p *PortfolioS) GetWeights() *ma.Matrix { return p.weights }
func (p *PortfolioS) SetWeights(w *ma.Matrix) { p.weights = w }
func (p *PortfolioS) GetAssets() *ma.Matrix { return p.assets }
func (p *PortfolioS) SetAssets(a *ma.Matrix) { p.assets = a }
func (p *PortfolioS) GetTickers() *ma.Matrix { return p.tickers }
func (p *PortfolioS) GetRiskReturnRatio() float64 { return p.risk / p.ret }

// temporary objects
func (p *PortfolioS) GetTempMeanReturns() *ma.Matrix { return p.rx }
func (p *PortfolioS) GetTempCovariances() *ma.Matrix { return p.cov }
func (p *PortfolioS) GetTempA() *ma.Matrix { return p.a }
func (p *PortfolioS) GetTempAinv() *ma.Matrix { return p.a_inv }
func (p *PortfolioS) GetTempB() *ma.Matrix { return p.b }
func (p *PortfolioS) GetTempX() *ma.Matrix { return p.x }
// 
// интерфес для метода наискорейшего спуска
//
// фукнция миниму которой мы ищем в методах многомерной оптимизации
func (p *PortfolioS) Function() float64 {
	var workers int
	var wg sync.WaitGroup
	//var Riski float64

	Reti := p.GetReturn()
	// получаем временые таблицы
	COV := p.GetTempCovariances()
	Ainv := p.GetTempAinv()
	B := p.GetTempB()
	X := p.GetTempX()
	// и список долей 
	W := p.GetWeights()

	//
	// Risk portrfolio variance
	// Ret portfolio return
	// Risk = f(Ret)
	//
	f := func (Reti float64) (Riski float64) {
		//
		// create B
		//
		for i := 0; i < Ainv.Rows()-2; i++ {
			B.Setij(i, 0, float64(0))
		}
		B.Setij(Ainv.Rows()-2, 0, float64(1))
		B.Setij(Ainv.Rows()-1, 0, float64(Reti))

		//
		// find weights as X = Ainv * B
		//
		ma.Mult(X, Ainv, B)

		//
		// calculate Risk
		//

		// get wights from X that we've found
		for i := 0; i < W.Rows(); i++ {
			value,_ := X.Getij(i, 0)
			W.Setij(i, 0, value)
		}

		// get sums for all weights
		// Sumi = Wi * Wj * COVj
		workers = W.Rows()
		wg.Add(workers)
		Tmp,_ := ma.New(W.Rows(), 1)
		for i := 0; i < workers; i++ {
			go func (a, b, c *ma.Matrix, cursor int, wg *sync.WaitGroup) {
				defer wg.Done()
				var sum float64 = 0.0
				// for each COV
				for j := 0; j < b.Cols(); j++ {
					value_wcursor,_ := c.Getij(cursor, 0)
					value_wj,_ := c.Getij(j, 0)
					value_cov,_ := b.Getij(cursor, j)
					sum += value_wcursor * value_wj * value_cov
					a.Setij(cursor, 0, sum)
				}
			}(Tmp, COV, W, i, &wg)
		} // eof for
		wg.Wait()

		//calc Risk as sum of weight sums
		Riski = ma.Sum(Tmp, 0, ma.COLVECTOR)

		//ma.Mprintf("\n##W", "%f", W)
		//ma.Mprintf("\n##TMP", "%f", Tmp)

		return
	}

	return f(Reti)
}
// вектор аргументов функции
func (p *PortfolioS) GetX() []float64 { return []float64{p.ret} }
func (p *PortfolioS) SetX(args []float64) { p.ret = args[0] }

// значение фукнции
func (p *PortfolioS) GetY() float64 { return p.risk }
func (p *PortfolioS) SetY(y float64) { p.risk = y }

// градиет функции - вектор частных производных
func (p *PortfolioS) Grad() []float64 {
	// get current args state
	args := p.GetX()
	// left and right deviations
	left,_ := New(p.GetAssets().Rows(), p.GetAssets().Cols())
	right,_ := New(p.GetAssets().Rows(), p.GetAssets().Cols())
	left.Assign(p)
	right.Assign(p)
	// numerical differentiation
	//
	// NOTE: delta ought to be small enough but you should remember 
	//       that too small value will drive to reducing accuracy
	//
	// df/dxi = f(x1,x2,...,xi+/\xi,...xn) - f(x1,x2,...xi-/\xi,...xn) / (2 * /\xi) 
	//
	delta := 0.05
	// vector of approx gradient of function f(x) values
	gradient := make([]float64, len(args))

	for i := 0; i < len(gradient); i++ {
		// differences
		plus := make([]float64,len(args))
		minus := make([]float64, len(args))
		copy(plus, args)
		copy(minus, args)
		// make a diffrence 
		plus[i] += delta
		minus[i] -= delta
		left.SetX(plus)
		right.SetX(minus)
		// calc aprox derivative
		gradient[i] = ( left.Function() - right.Function() ) / (2.0 * delta)
	}
	return gradient
}

// задать фукнцию оценки g(x) для метода одномерной оптимизации
func (p *PortfolioS) GetG() (func (data ma.ODO, lambda float64) float64) {
	g := func (data ma.ODO, lambda float64) float64 {
		args := data.GetX()
		estimation,_ := New(p.GetAssets().Rows(), p.GetAssets().Cols())
		estimation.Assign(p)
		// vector of approx gradient of function f(x) values
		gradient := data.Grad()
		for i,df_dxi := range gradient {
			diff := make([]float64, len(args))
			copy(diff, args)
			diff[i] = args[i] - lambda * df_dxi
			estimation.SetX(diff)
		}
		return estimation.Function()
	}
	return g
}

// отрезок для поиска минимума методом одной оптимизации a,b 
func (p *PortfolioS) GetInterval() (float64,float64) {
	a, b := -10000.0, 100000.0
	return a,b
}
// максимальное количесво приближений 
//n int - для одномерной оптимизации
func (p *PortfolioS) GetMaxIterODO() (n int) {
	n = 100
	return
}
//m int - многомерной оптимизации
func (p *PortfolioS) GetMaxIterMDO() (m int) {
	m = 1000
	return
}
// требуемая точность
//epsilon1 float64 для метода одномерной оптимизации
func (p *PortfolioS) GetEpsilon1() float64 { return 0.01 }
//epsilon2 float64 для метода многомерной оптимизации
func (p *PortfolioS) GetEpsilon2() float64 { return 0.01 }
// Аргумент минимизации (для наискорейшего спуска argmin по lambda)
// phi(x) метод одномерной оптимизаци
func (p *PortfolioS) GetF() (func(data ma.ODO) float64) {
	return ma.Dichotomia
	//return ma.Goldensection
}




//
// Интерфейс sort.Sort() для популяции
//
type PopulationS struct {
	individuals []*PortfolioS
	assets *ma.Matrix
}
func NewPop(quantity int) *PopulationS {
	pop := new(PopulationS)
	pop.individuals = make([]*PortfolioS, quantity)
	pop.assets = nil
	return pop
}
// PopulaionS sort.Sort() interface methods
// Swap methods of the embedded type value.
// exported sort.Interface Len() int
func (p *PopulationS) Len() int { return len(p.individuals) }
// Swap methods of the embedded type value.
// exported sort.Interface Swap(i,j int)
func (p *PopulationS) Swap(i, j int) { p.individuals[i], p.individuals[j] = p.individuals[j], p.individuals[i] }
// Methods implements sort.Interface by providing 
// Less and using the Len and Swap methods of the embedded type value.
// sort descent by risk return ratio
//type ByRiskReturnRatio struct{ PopulationS }
//func (p *ByRiskReturnRatio) Less(i, j int) bool {
//	return p.PopulationS[i].GetRiskReturnRatio() < p.PopulationS[j].GetRiskReturnRatio()
//}
func (p *PopulationS) Less(i, j int) (res bool) {
	return p.individuals[i].GetReturn() < p.individuals[j].GetReturn()
}

// population interface
// исходная популяция или начальное приближение
// отец-основатель сущ founding father
func (p *PopulationS) FoundingFathers() {
	var workers int
	var wg sync.WaitGroup
	//
	// загрузить данные о доступных инструментах
	//
	fmt.Println("FoundingFathers()")
	//FIXME: open file for reading
	fi, err := os.Open("./t_assets_20.mat")
	if err != nil {
		fmt.Println("Error openning file:", err)
		os.Exit(1)
	}
	defer fi.Close()

	// read file line by line
	matrices, err := ma.Fnew(fi)
	if err != nil {
		fmt.Println("Error read file:", err)
		os.Exit(1)
	}

	r := matrices[0]
	ma.Mprintf("\n#", "%.2f ", r)
	Assets,_ := ma.New(r.Rows(),r.Cols())
	ma.Assign(Assets,r)

	//
	// получить перивичные параметры
	//
	//количество инструментов в портфеле
	quantity := 5
	// портфель состоит из 5 активов всего доступных активов 20

	//
	// сформировать поколоение начального приближения
	//
	workers = p.Len()
	wg.Add(workers)
	for i := 0; i < p.Len(); i++ {
		// запустить воркера для инициализации и формализации особи 
		go func (p *PopulationS, a *ma.Matrix, index, quantity int, wg *sync.WaitGroup) {
			defer wg.Done()
			portfolio,_ := New(a.Rows(), quantity)
			portfolio.Make(a)
			portfolio.Formalize()
			p.individuals[index] = portfolio
		}(p, Assets, i, quantity, &wg)
	} // eof for
	wg.Wait()
	//
	// запомнить исходную таблицу
	//
	p.assets = Assets
	fmt.Printf(" Assets: %d\n", Assets.Cols())
}
// вычисление функции пригодности полученных решений
func (p *PopulationS) EvaluateFitness() {
	//fmt.Println("EvaluateFitness()")
	var workers int
	var wg sync.WaitGroup
	// для каждого портфеля
	// найти оптимальный портфель 
	// по условию оптимальности:
	//   максимальная доходность
	//   при минимальном риске
	//   при этом требуется неотрицательные веса
	//   и равномерное распределение по инструментам(опция)
	workers = p.Len()
	wg.Add(workers)
	for i := 0; i < p.Len(); i++ {
		// запустить воркера для определения пригодности особи
		go func (p *PopulationS, index int, wg *sync.WaitGroup) {
			defer wg.Done()
			// получить оптимальную конфигурацию портфеля по соотноешнию Risk/Return
			p.optimize(index)
			// Penalties
			// тут надо сделать систему штрафов для отбора действительно пригодных решений
			// оценить полученный оптимальный портфель
			p.estimate(index)
		}(p, i, &wg)
	} // eof for
	wg.Wait()
}
// отбор наиболее пригодных решений
func (p *PopulationS) SelectFittestSurvivors() {
	//fmt.Println("SelectFittestSurvivors()")
	// произвести отбор, здесь просто сортировка
	// лучшие будут сосредоточены с одного из концов списка
	sort.Sort(p)

	// DEBUG
	//for i := 0; i < p.survivors(); i++ {
	//	portfolio := p.individuals[i]
	//	fmt.Printf("\n DOOMED[%d]:Return: %.3f Risk: %.3f Risk/Return: %.3f\n", i, portfolio.GetReturn(), portfolio.GetRisk(), portfolio.GetRiskReturnRatio())
	//}
	//fmt.Println()

	// DEBUG
	for i := p.survivors(); i < p.Len(); i++ {
		portfolio := p.individuals[i]
		fmt.Printf("\n SURV[%d]:Return: %.3f Risk: %.3f Risk/Return: %.3f\n", i, portfolio.GetReturn(), portfolio.GetRisk(), portfolio.GetRiskReturnRatio())
	}
	fmt.Println()

	// show fittest
	portfolio := p.individuals[p.Len()-1]
	ma.Mprintf("  Portfolio:", "%0.f", portfolio.GetTickers())
	ma.Mprintf("  R:", "%f", portfolio.GetTempMeanReturns())
	ma.Mprintf("  W:", "%f", portfolio.GetWeights())
	fmt.Printf("\n FITTEST[%d]:Return: %.3f Risk: %.3f Risk/Return: %.3f\n", p.Len()-1, portfolio.GetReturn(), portfolio.GetRisk(), portfolio.GetRiskReturnRatio())
	fmt.Println()
}
// Выполнить операцию рекомбинации. Сформировать новое множество 
// решений из наиболее пригодных при помощи генетических операторов
func (p *PopulationS) Recombine() {
	//fmt.Println("Recombine()")
	var workers int
	var wg sync.WaitGroup
	// Оператор рекомбинации - языковая конструкция, которая определяет, 
	// как новая генерация хромосом будет построена из родителей и потомков.
	// Важным понятием при реализации генетических операторов является вероятность,
	// которая определяет, какой процент общего числа генов в популяции изменяется каждой генерации.
	//TODO: Для оптимизационных задач обычно принимают следующие вероятности операторов:
	//   вероятность кроссинговера обычно принимают в пределах (0,6¸0,99);
	//   мутации – 0,6;
	//   инверсии – (0,1¸0,5);
	//   транслокации - (0,1¸0,5);
	//   транспозиции - (0,1¸0,5);
	//   сегрегации - (0,6¸0,99);
	//   удаления - (0,6¸0,99);
	//   вставки - (0,6¸0,99).
	// Мои операторы
	//   репарации - (0,3?)
	//   миграции - (0.1?)
	workers = p.survivors()
	wg.Add(workers)
	for i := 0; i < p.survivors(); i++ {
		// запустить воркера для создания нового индивида
		go func (p *PopulationS, index int, wg *sync.WaitGroup) {
			defer wg.Done()
			// разыграть шанс для того или оного оператора
			chance := ma.RndBetweenU(0,100)
			//fmt.Printf("GA operation chance: %d operation: ", chance)
			if 85 < chance && chance <= 100 {
			//80..100
				p.repair(index)
			} else if 75 < chance && chance <= 85 {
			//70..80
				p.migrate(index)
			} else if 60 < chance && chance <= 75 {
			//60..70
				p.inverse(index)
			} else if 40 < chance && chance <= 60 {
			//40..60
				p.mutate(index)
			} else {
			//0..40
				p.crossing_over(index)
			}
		}(p, i, &wg)
	}
	wg.Wait()
}
// оптимизация портфеля
func (p *PopulationS) optimize(index int) {
	portfolio := p.individuals[index]
	// задать исходное приближение
	Rp0 := 0.5
	portfolio.SetReturn(Rp0)
	ma.SteepestDescent(portfolio)
	//fmt.Printf("\n[%d]:Return: %.3f Risk: %.3f Risk/Return: %.3f\n", index, portfolio.GetReturn(), portfolio.GetRisk(), portfolio.GetRiskReturnRatio())
	//ma.Mprintf("  portfolio: ", "%.0f", portfolio.GetTickers())
}
// оценка полученного оптимального портфеля
func (p *PopulationS) estimate(index int) {
	portfolio := p.individuals[index]
	// Штрафная функция
	// назначить штраф за отрицательный вес в портфеле
	minWeight,_ := ma.Min(portfolio.GetWeights(),0,ma.COLVECTOR)
	if minWeight < 0 {
		portfolio.SetReturn(0.01)
	}
	// Штрафная функция
	// назначить штраф за повторяющиеся инструменты(дубли)
	// Простой однопроходный алгоритм поиска одинаковых элементов списка(массива)
	// основанный на использовании хэш-таблицы для определения(и, возможно,хранения) дублей
	test_duplicate_items := func (list *ma.Matrix) (multiple bool) {
		// хэш таблица для определения(и,возможно, хранения) дублей
		hash := make(map[float64]bool)
		multiple = false
		//stop := false
		// проход по списку(массиву) с целью выявленяи дублей
		for i := 0; i < list.Rows(); i++ {
			// получить элемент списка
			item,_ := list.Getij(i,0)
			// проверить существование записи с ключом элемента в хэш-таблице
			// если ключ такой есть - элемнт уже встречалься, то останов
			if _,ok := hash[item]; ok {
				multiple = true
				//stop := true
				break
			}
			// внести в хэш-таблицу запись с ключом в виде элемента списка
			hash[item] = true
		}
		return
	}
	if test_duplicate_items(portfolio.GetTickers()) {
		portfolio.SetReturn(0.01)
	}
}
// размер выбырки наиболее пригодных
func (p *PopulationS) survivors() (survivors int) {
	survivors = p.Len() - int(p.Len() / 4) - 1
	if survivors <= 0 {
		survivors = 1
	}
	return
}
// получить родителя из наиболее пригодных особоей
func (p *PopulationS) getParent() *PortfolioS {
	parent := ma.RndBetweenU(p.survivors(), p.Len())
	return p.individuals[parent]
}
// Генетические операторы(свободные вариации на тему)
// [Оператор миграции]
// это не совсем оператор, просто перенос генов без измениний, миграция особи из поколения в поколение
func (p *PopulationS) migrate(index int) {
	//fmt.Println("migrate")
	//fmt.Println("-----------START---------------")
	// Выбрать одного родителея для создания новой особи из наиболее пригодных
	portfolio_A := p.getParent()
	// сформировать ребнка
	portfolio_B,_ := New(portfolio_A.Rows(), portfolio_A.Cols())
	// Передать генотип выбраной особи новой особи без внесения измений  
	portfolio_B.Assign(portfolio_A)
	// привести портфель в надлежащий вид для решения задачо оптимизации
	portfolio_B.Formalize()
	//ma.Tprintf("  Parent:", "%.2f", portfolio_A.GetTickers(), portfolio_A.GetTempMeanReturns())
	//ma.Tprintf("  Child:", "%.2f", portfolio_B.GetTickers(), portfolio_B.GetTempMeanReturns())
	// поместить новую особь в популяцию нового поколения
	p.individuals[index] = portfolio_B
	//fmt.Println("-------------END---------------")
}


// [Универсальный оператор кроссинговера]
func (p *PopulationS) crossing_over(index int) {
	//fmt.Println("crossing_over")
	//fmt.Println("-----------START---------------")
	// Выбрать родителей A и B для создания новой особи из наиболее пригодных
	portfolio_A := p.getParent()
	portfolio_B := p.getParent()
	// сформировать ребнка
	portfolio_C,_ := New(portfolio_A.Rows(), portfolio_A.Cols())
	// Передать генотип одного из радителей для инициализации(здесь: родитель А)
	portfolio_C.Assign(portfolio_A)
	// выбрать в каких долях гены родителей будут присуствовать в потомке
	// доля не может быть меньше 20% и не больше 80%
	racio := ma.RndBetweenU(20,80)
	//fmt.Printf("PrantA portion: %d%% ParentB portion: %d%%\n", racio, 100-racio)
	// для каждого ген в хромосоме разыграть шанс получения гена первого или вртого родителя
	for gene := 0; gene < portfolio_C.GetTickers().Rows(); gene++ {
		chance := ma.RndBetweenU(0,100)
		//fmt.Printf("[%d] Chance: %d ", i, chance)
		if chance <= racio {
		// parent A gene
		//	fmt.Printf(" A")
			portfolio_C.Replace(portfolio_A, gene)
		} else {
		// parent B gene
		//	fmt.Printf(" B")
			portfolio_C.Replace(portfolio_B, gene)
		}
	}
	//fmt.Println()
	// привести портфель в надлежащий вид для решения задачо оптимизации
	portfolio_C.Formalize()
	//ma.Tprintf("  Parent_A:", "%.2f", portfolio_A.GetTickers(), portfolio_A.GetTempMeanReturns())
	//ma.Tprintf("  Parent_B:", "%.2f", portfolio_B.GetTickers(), portfolio_B.GetTempMeanReturns())
	//ma.Tprintf("  Child:", "%.2f", portfolio_C.GetTickers(), portfolio_C.GetTempMeanReturns())
	// поместить новую особь в популяцию нового поколения
	p.individuals[index] = portfolio_C
	//fmt.Println("-------------END---------------")
}

// [Оператор мутации]
func (p *PopulationS) mutate(index int) {
	//fmt.Println("mutate")
	//fmt.Println("-----------START---------------")
	// Выбрать одного родителея для создания новой особи из наиболее пригодных
	portfolio_A := p.getParent()
	// сформировать ребнка
	portfolio_B,_ := New(portfolio_A.Rows(), portfolio_A.Cols())
	// Передать генотип выбраной особи новой особи без внесения измений  
	portfolio_B.Assign(portfolio_A)
	// определяется позиция мутирующего гена
	gene := ma.RndBetweenU(0,portfolio_A.GetTickers().Rows()) //от 0 до Количество генов в хромосоме
	// определяется новый параметр гена
	ticker := ma.RndBetweenU(0,p.assets.Cols()) // дипазон достуных значений(зависит от реализации)
	//fmt.Printf("Gene: %d  tikcer: %d\n", gene, ticker)
	portfolio_B.Take(p.assets, ticker, gene)
	// привести портфель в надлежащий вид для решения задачо оптимизации
	portfolio_B.Formalize()
	//ma.Tprintf("  Parent:", "%.2f", portfolio_A.GetTickers(), portfolio_A.GetTempMeanReturns())
	//ma.Tprintf("  Child:", "%.2f", portfolio_B.GetTickers(), portfolio_B.GetTempMeanReturns())
	// поместить новую особь в популяцию нового поколения
	p.individuals[index] = portfolio_B
	//fmt.Println("-------------END---------------")
}

// [Оператор инверсии]
func (p *PopulationS) inverse(index int) {
	//fmt.Println("inverse")
	//fmt.Println("-----------START---------------")
	// Выбрать одного родителея для создания новой особи из наиболее пригодных
	portfolio_A := p.getParent()
	// сформировать ребнка и придать емукачества первого родителя
	portfolio_B,_ := New(portfolio_A.Rows(), portfolio_A.Cols())
	// Передать генотип
	portfolio_B.Assign(portfolio_A)
	//получить инвертируемый сегмент хромосомы
	y1 := float64(ma.RndBetweenU(0,portfolio_A.Cols()))
	y2 := float64(ma.RndBetweenU(0,portfolio_A.Cols()))
	//определить левую и правую границы
	a := int(math.Min(y1,y2))
	b := int(math.Max(y1,y2))
	//произвести инвертирование генов на указанном диапазоне
	//fmt.Printf(" a: %d  b: %d\n", a, b)
	for gene := a; gene < b; gene++ {
		// здесь эта операция чень похожа на мутацию но несколько иная по реализации
		// и производит изменения на участке a,b
		ticker,_ := portfolio_A.GetTickers().Getij(gene, 0)
		ticker_inv := p.assets.Cols() - int(ticker)
		//fmt.Printf("[%d] value: %.0f  -> value_inv: %d\n", gene, ticker, ticker_inv)
		portfolio_B.Take(p.assets, ticker_inv, gene)
	}
	// привести портфель в надлежащий вид для решения задачо оптимизации
	portfolio_B.Formalize()
	//ma.Tprintf("  Parent:", "%.2f", portfolio_A.GetTickers(), portfolio_A.GetTempMeanReturns())
	//ma.Tprintf("  Child:", "%.2f", portfolio_B.GetTickers(), portfolio_B.GetTempMeanReturns())
	// поместить новую особь в популяцию нового поколения
	p.individuals[index] = portfolio_B
	//fmt.Println("-------------END---------------")
}

// [Оператор репарации]
func (p *PopulationS) repair(index int) {
	//fmt.Println("repair")
	//fmt.Println("-----------START---------------")

	// здесь эта функция реализована путём внедрения в новое поколение "середняков" из прошлого поколения
	// "средняя особь" выбирается как значение медианы

	// Выбрать одного родителея для создания новой особи как значение медианы популяции
	// медианна будет соответствовать центральному значению ряда
	// (без учёта чётности количества елементов, здесь это не так важно),
	// номер которого можно определить по формуле:
	// IndexMe = N/2, где
	// IndexMe - номер значения, воответствуюцего медиане,
	// N - количество знаений в совокупности данных.	
	indexMe := int(p.Len() / 2)
	//fmt.Printf("pop len: %d  indexMe: %d\n", p.Len(), indexMe)
	portfolio_A := p.individuals[indexMe]
	// сформировать ребнка
	portfolio_B,_ := New(portfolio_A.Rows(), portfolio_A.Cols())
	// Передать генотип выбраной особи новой особи без внесения измений  
	portfolio_B.Assign(portfolio_A)
	// привести портфель в надлежащий вид для решения задачо оптимизации
	portfolio_B.Formalize()
	//ma.Tprintf("  Parent:", "%.2f", portfolio_A.GetTickers(), portfolio_A.GetTempMeanReturns())
	//ma.Tprintf("  Child:", "%.2f", portfolio_B.GetTickers(), portfolio_B.GetTempMeanReturns())
	// поместить новую особь в популяцию нового поколения
	p.individuals[index] = portfolio_B

	//fmt.Println("-------------END---------------")

	// Немножго теории:
	// Формула медианы для дискретных данных чем-то напоминает формулу моды. 
	// А именно тем, что формулы как таковой нет. Медианное значение выбирают 
	// из имеющихся данных и только, если это невозможно, проводят несложный расчет.

	// Первым делом данные ранжируют (сортируют по убыванию). 
	// Далее есть два варианта. Если количество значений нечетно, 
	// то медианна будет соответствовать центральному значению ряда, 
	// номер которого можно определить по формуле:
	// IndexMe = (N+1)/2, где
	// IndexMe - номер значения, воответствуюцего медиане,
	// N - количество знаений в совокупности данных.	
	// Тогда медиана будет обозначаться, как
	// Me = (XN+1)/2
	// Определение медианы по центральному значению
	// Это первый вариант, когда в данных есть одно центральное значение. Второй вариант наступает тогда, когда количество данных четно, то есть вместо одного есть два центральных значения. Выход прост: берется средняя арифметическая из двух центральных значений:
	// Me = (XN/2 + XN/2+1)/2
	// Определение медианы при четном количестве данных
	// Так происходит поиск или расчет в дискретных данных.

	/*

	[https://laservirta.ru/%D1%84%D0%BE%D1%80%D0%BC%D1%83%D0%BB%D0%B0-%D0%BC%D0%BE%D0%B4%D1%8B-%D0%B8-%D0%BC%D0%B5%D0%B4%D0%B8%D0%B0%D0%BD%D1%8B-%D0%B2-%D1%81%D1%82%D0%B0%D1%82%D0%B8%D1%81%D1%82%D0%B8%D0%BA%D0%B5/]

	Медиану следует применять в качестве средней величины в тех случаях, 
	где нет достаточной уверенности в однородности изучаемой совокупности. 
	На медиану влияют не столько сами значения, сколько число случаев на том 
	или ином уровне. Следует также отметить, что медиана всегда конкретна 
	(при большом числе наблюдений или в случае нечетного числа членов совокупности), 
	т.к. под Ме подразумевается некоторый действительный реальный элемент совокупности, 
	тогда как арифметическая средняя часто принимает такое значение, которое не может 
	принимать ни одна из единиц совокупности.
	

	Главное свойство Ме в том, что сумма абсолютных отклонений значений признака от 
	медианы меньше, чем от любой другой величины:  
	sum(i=1,n)|xi - Me| = min
	. Это свойство Ме может быть использовано, 
	например, при определении места строительства общественных зданий, т.к. Ме определяет точку, 
	дающую наименьшее расстояние, допустим, детских садов от местожительства родителей, 
	жителей населенного пункта от кинотеатра, при проектировке трамвайных, троллейбусных остановок и т.д.
	В системе структурных показателей в качестве показателей особенностей формы распределения выступают варианты, 
	занимающие определенное место в ранжированном вариационном ряду (каждое четвертое, пятое, десятое, двадцать пятое и т.д.). 
	Аналогично с нахождением медианы в вариационных рядах можно отыскать значение признака у любой по порядку единицы ранжированного ряда.

	

	посмотреть квартили и децили
	*/
}

// условие останова
func (p PopulationS) Done(generation int) bool {
	// если работа только начата, то проверять условия останова не следует
	if generation < 5 {
		return false
	}
	//TODO: эпсилон надо бы задавать из вне
	epsilon := 0.001
	squareMean := 0.0
	for i := p.Len()-1; i > p.survivors()-1; i-- {
		portfolio_A := p.individuals[i]
		portfolio_B := p.individuals[i-1]
		squareMean += (portfolio_A.GetReturn() - portfolio_B.GetReturn())*(portfolio_A.GetReturn() - portfolio_B.GetReturn())
	}
	fmt.Printf(" Norm: %f\n", squareMean)
	if squareMean <= epsilon {
		return true
	} else {
		return false
	}
}

//
// main driver
//
func main() {


	// shake the generator!
	ma.SRnd64(time.Now().Unix())

	// 
	// Genetic Algorithm
	//

	// создать популяцию численности N
	populationSize := 100
	pop := NewPop(populationSize)

	//ma.GeneticAlgorithm(pop)

	// сформировать исходную популяцию или начальное приближение
	pop.FoundingFathers()

	// максимальное число итераций
	var niterations = 100
	for i := 0; i < niterations && !pop.Done(i); i++ {
		fmt.Printf("  GENERATION[%d]\n", i)

		// Estimate population fitnesses
		pop.EvaluateFitness()

		// There will be selection of the fittest
		pop.SelectFittestSurvivors()

		// Spawn new generation from the most fittest
		pop.Recombine()
		//fmt.Printf("\n")
	}
	//pop.Result()


} // eof main


