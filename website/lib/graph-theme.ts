/** Presentation-only colors; never encode a shared biological score. */
export const graphTheme = {
  background: '#f4f8fc',
  ink: '#294760',
  muted: '#668097',
  faintEdge: '#bacbdb66',
  selection: '#345fac',
  colors: {
    candidate: '#269c90', anchor: '#527bcc', both: '#a477ba',
    protein: '#3d9db5', term: '#9e84bd', seed: '#bc914b',
    experimental: '#38a496', predicted: '#6a8ccc', isoform: '#aa86c5',
    annotation: '#9aa8c4', omnipath: '#62a9b6', homology: '#bc914b',
  } as Record<string, string>,
  density: ['#eaf2f7', '#c8e4e6', '#89c9c5', '#47a7a7', '#267f8e', '#21546f'],
};
