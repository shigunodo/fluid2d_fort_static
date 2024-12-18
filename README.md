本プログラムは、Qitta の記事

+ [流体計算の実行速度比較: Fortran, C++, Rust, Python, Julia](https://qiita.com/shigunodo/items/d693dc03323f9a205bb9)

での実行速度比較用途で作成した。同記事の「Fortran/静的」の場合のコード/設定である。 

## 計算内容

詳しくは上で示した記事に書かれているが、2 次元圧縮性オイラー方程式の初期値問題（時間発展問題）を解くプログラムである。一般の構造格子上で風上差分法を用いて解かれる。

### 記事で説明されていないこと

+ 計算スキームは記事で用いた Roe型 FDS 法 & MP5 法以外にも実装されているが、デフォルトでは記事の設定になっている。スキームの変更はファイル ``proceed.f90`` 内のサブルーチン ``RHS`` 内で、変数 ``con_sch``, ``ma_sch`` の変更によって行う。
+ 計算領域の境界には、境界条件の反映のために、ダミーのグリッドが NB 個存在する。これはグリッド数に含まれている。つまり、i 方向に NI 個、j 方向に NJ 個の構造格子上で解く場合、真に方程式の右辺が評価されるのは (NI - 2 \* NB) \* (NJ - 2 \* NB) 個のグリッドに限られる。

## 動かし方

### 座標と初期条件を与える

プログラムを動かすためには、グリッド座標の一覧と初期条件を外部ファイルとして与える必要がある。

i 方向にグリッド数 NI、j 方向にグリッド数 NJ の構造格子上で解く場合、各量は NI \* NJ 配列で与えられる。各配列は Column-major の順序で 1 次元的にリシェイプされ、改行によって縦に並べる書式で外部ファイル内に記述される必要がある。

外部ファイルが置かれるべきディレクトリ (ワーキングディレクトリ) のパス指定はファイル ``parameters_io.f90`` 内で宣言された変数 ``dir`` によって行う。

+ **グリッド座標**: ワーキングディレクトリ下に ``coordinate.dat`` というファイル名で置かれるべき。各グリッドの x 座標を表す配列の後に、続けて y 座標を表す配列が記述されているべき。
+ **初期条件**: ワーキングディレクトリ下に、``b0000000.dat`` という名前で置かれているべき。流体の密度 (rho)、x 流速 (u)、y 流速 (v)、単位体積当たりの全エネルギー (e) の順で記述されているべき。

### コンパイルする

例えば、``gfortran`` でコンパイルする場合、``src/`` ディレクトリ下で

```bash
gfortran main.f90
```

とすれば、実行ファイルが生成される。

実行すると、``parameters_io.f90`` 内で宣言された変数 ``Tmax`` だけ時間積分が行われる。変数 ``Nout`` に指定した数だけ、等物理時間間隔で計算結果がワーキングディレクトリ下に出力される。ファイル出力時には、標準出力にステータスが追加表示される。ステータスと計算設定はファイル ``settings.dat`` にも出力される。

### 座標と初期条件の生成

``data/`` ディレクトリに、記事の KH 不安定計算用の座標・初期条件を置いておいた。

また、同ディレクトリには ``fluid.py`` という Python スクリプトファイルも置かれている。これは、KH 不安定用の格子・初期条件を行ったり、計算結果をプロットするためのものである。

```Python
import fluid
Ni = 408 # i 方向グリッド数
Nj = 408 # j 方向グリッド数
Nb = 4 # 境界条件用グリッド数
gamma = 1.4 # 比熱比
a = fluid.Fluid2d(Ni,Nj,Nb,gamma)
```

などとすれば、内部に KH 不安定用の座標と初期条件を持ったオブジェクト ``a`` が生成される。

```Python
a.output_coordinate_fort("座標ファイルパス")
a.output_basic_fort("初期条件ファイルパス")
```

とすれば、それをファイル出力できる。

計算結果をプロットしたい場合は、

```Python
a.input_basic_fort("ファイルパス")
```

によって結果を読み込んでから、

```Python
a.plot_rho("画像ファイルパス.png", time) # time はタイトルに表示する物理時間 (実数)
```

とすれば、密度をカラーマップに図示した png 画像が保存される。