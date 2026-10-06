import {lazy, Suspense, type ComponentType} from 'react';

/** Browser-only lazy component, including React 19 ref-as-prop forwarding. */
export default function dynamic<P extends object>(
  loader: () => Promise<{default: ComponentType<P>}>,
  options: {ssr?: boolean; loading: ComponentType},
) {
  const Component = lazy(loader);
  const Loading = options.loading;
  return function LazyClient(props: P) {
    return <Suspense fallback={<Loading />}><Component {...props} /></Suspense>;
  };
}
