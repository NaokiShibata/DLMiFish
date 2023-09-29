# DLMiFish (Beta version)
NCBIから真骨類 (Teleostei, taxid=[32443](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=32443&lvl=3&lin=f&keep=1&srchmode=1&unlock))、軟骨魚類 (Chondrichthyes, taxid=[7777](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=7777&lvl=3&lin=f&keep=1&srchmode=1&unlock))、円口類(Cyclostomata, taxid=[1476529](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1476529&lvl=3&lin=f&keep=1&srchmode=1&unlock))の12s rRNA領域の配列データをダウンロードするパッケージです。  
`Biopython`に搭載されている`Entrez`の各関数を使って、上記3分類群のGenbankの情報を取得した後に12s rRNA領域の配列情報を抽出します。また、Optionでは、`Cutadapt`を使ってプライマー配列の一部削除(データ軽量化)や`seqkit`による配列データの配列長を視覚化もできたりします。  
1度ダウンロードした.gbファイルはGI_Folderに保存されます。ダウンロードする.gbファイルはNCBIでの問い合わせによって得られたACCの情報と、これまでにダウンロードした.gbファイル名に使用されているACCの情報を比較し、今回とこれまでの差分をとることでダウンロードする.gbファイルを決めています。

windows10 WSLとUbuntu 20.04 LTSでの動作を確認しています。  

## 動作環境や環境整備に必要なもの

- python3.6 以上 (f string使用のため)  
以下は最新版 or pythonのバージョンにあったもの  
  - pandas
  - biopython
  - tqdm
  - joblib
  - cutadapt
- conda : 最新版 or pythonのバージョンにあったもの
- [mamba](https://github.com/mamba-org/mamba) : 最新版 or pythonのバージョンにあったもの
- seqkit : version2.0以上

## 01 仮想環境の作成
pythonでは、メイン環境とは他に仮想的な実行環境を作成できます。  

主な目的は以下の通り
- システム全体で使うpython環境に影響を与えずモジュールの追加や入れ替えをしたい
- 異なるversionのpythonを使い分けたりしたい
- 失敗した場合に、簡単に削除できる  

pythonの仮想環境はPython3 の標準ライブラリである`venv`や`virtualenv`などがあります。  
minicondaをベースとした`Conda`での仮想環境作成法は以下のブログ記事にて紹介していますのでご参照ください。  
- [pythonの環境構築とその他もろもろについての覚え書き](https://edna-blog.com/technique/python_env/) ＠はじめての環境DNA  

## 02モジュールをインストールする
まず、パッケージマネージャーである`mamba`をインストールします。mambaはマルチスレッドでリポジトリのデータやパッケージファイルを並列ダウンロードしたり、依存関係を管理してくれたりします。

- [mamba](https://github.com/mamba-org/mamba) : Github repository  

  ### Installing based on conda
  ```bash
  conda install mamba -n base -c conda-forge
  ```
その他、pythonモジュールのインストール
- [biopython](https://github.com/biopython/biopython) : Github repository  
- [pandas](https://github.com/pandas-dev/pandas) : Github repository  
- [tqdm](https://github.com/tqdm/tqdm) : Github repository  

  ### Install
  ```bash
  mamba install bioconda bipython cutadapt -y
  mamba install -c anaconda pandas -y
  mamba install -c conda-forge tqdm joblib -y
  ```
次に、fastqやfastaの操作ツールである`seqkit`をインストールします。 

- [seqkit](https://github.com/shenwei356/seqkit) : Github repository  
  ### Install
  __Windows 10でcondaを使わない場合__
  
  1. Homeデレクトリに移動
  ```bash
  cd
  ```
  
  2. 最新版をダウンロード(要version確認)
  ```bash
  # 2022/4/21
  wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_arm64.tar.gz
  ```
  
  3. 解凍
  ```
  tar -zxvf *arm64.tar.gz
  ```
  
  4. プログラムの置き場に移動
  ```
  sudo mv seqkit ~/usr/local/bin/
  ```
  
  __`conda`が使える方は以下で一発インストール__
  ```bash
  conda install -c bioconda seqkit -y
  ```
  
  __`mamba`でもよい__
  ```bash
  mamba install -c bioconda seqkit -y
  ```
  
  5. 確認
  ```bash
  seqkit version
  # seqkit v2.0.0
  ```
  
## 03 DLMiFishのダウンロード
デスクトップに配置することを想定

### using git
`git clone`を使用したパッケージ取得。デスクトップをカレントデレクトリとした後に以下を実行
```bash
git clone https://github.com/NaokiShibata/DLMiFish.git
```
### Download ZIP
Download ZIPボタンよりダウンロードすると`ZIP`ファイルがダウンロードされるので、デスクトップに置き、解凍します。
![image](https://user-images.githubusercontent.com/53568847/164349138-3227f1cd-3e16-45c8-a3e7-4868dfeb303c.png)

解凍したら、DLMiFish-masterというフォルダ名をDLMiFishに変更します。

## 04 使い方
### クイックスタート
`setting.txt`に以下の情報を追加して、`python3 DLMiFish.py`を実行すると、配列取得が開始されるはずです。  
注意点として、現状では1GB程度のファイルサイズになることと、環境にもよりますが、デフォルトの設定では3時間程度かかります(2022/4現在のデータ量)。

__APIs__  
NCBIのAPI keyを取得してください。Entrezの使用上は任意ですが、プログラムの作成上表記無しで動く設計にできていません。  
また、利用頻度が高い場合(3回/1秒以上)には必要になるので登録&取得しておきましょう。  

- [NCBI account](https://www.ncbi.nlm.nih.gov/account/)

アカウント登録とAPIの取得ができましたら、`seatting.txt`の下記部分に必要事項を記載します。
- `Email = e-mail adress` : 有効なe-mailアドレスを指定
- `Api_key =  API` : 取得したAPIを指定

注意点として、'='の後のemailとAPI keyは特に””や’’でくくる必要はありません。

その他、特筆すべき部分について記載します。

__取り出す情報数の指定__  
任意のワード検索によってHitした情報からいくつ取り出すかについて、`Count_in_query1(=retmax)`に指定可能です。

例 : Count_in_query1 = 10000 (最大数は10000)

__1回の要求あたりの情報のダウンロード数__  

例 : Count_in_query2 = 20 (おすすめは20程度)

#### オプションの設定
プライマー配列を認識して、プライマー配列とその外側を除去したり、配列長の閾値を設定したり出来ます。
オプションを有効化するには、`option = Yes`を指定してください。そうすることで以下のオプションが有効になります。

__プライマー配列より外側の配列の除去__  
`Cutadapt`による配列除去が実行されます。  
- `error = 0.1` : 配列一致度の許容率(プライマー配列長×指定値)
- `core = 24` : 使用するスレッド数(増やしてもあまり恩恵を感じれませんが、多めに指定しておいたら良いと思います)

__配列長によるフィルタリング__  
`seqkit`による配列長フィルタリングが実行されます。  
- `minilength = 120` : 可変領域の最小値をイメージして指定します。MiFishプライマーの可変領域の平均配列帳はca172bp、最小は140bp程度です。 

#### 実行

以下ですぐに実行が可能です。

```bash
python3 DLMiFish.py
```

※ 2023/9/30時点で*18041*本の配列がHitし、ダウンロードには50分かかりました。
