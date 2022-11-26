import { ReactChild } from 'react';
import { useOnboarding } from '../../hooks/onboarding';
import { MainLayoutRoute } from '../MainLayout';
import { TopBar } from '../TopBar';

export function OnboardingLayout({ children }: { children: ReactChild }) {
  const { state } = useOnboarding();

  return !state.alonePage ? (
    <>
      <TopBar progress={state.progress}></TopBar>
      <div className="flex-grow py-10 mx-4 h-full">{children}</div>
    </>
  ) : (
    <MainLayoutRoute widgets={false}>
      <div className="flex-grow pt-10 mx-4">{children}</div>
    </MainLayoutRoute>
  );
}
